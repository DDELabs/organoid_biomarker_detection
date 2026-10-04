"""Read a Bioconductor ExpressionSet from an .rda using the low-level rdata parser."""
import numpy as np, pandas as pd, rdata
from rdata.parser import RObjectType as T

def deref(o):
    while o is not None and o.info.type == T.REF:
        o = o.referenced_object
    return o

def sym(tag):
    tag = deref(tag)
    return tag.value.value.decode() if tag.info.type == T.SYM else str(tag.value)

def pairlist(o):
    out = {}
    while o is not None and o.info.type in (T.LIST, T.ATTRLIST if hasattr(T,'ATTRLIST') else T.LIST):
        car, cdr = o.value
        out[sym(o.tag)] = car
        o = cdr
    return out

def attrs(o):
    return pairlist(o.attributes) if o.attributes is not None else {}

def strings(o):
    o = deref(o)
    return [None if x.value is None else x.value.decode("latin1") for x in o.value]

def env_vars(env):
    env = deref(env)
    # ENV value: (enclos, frame, hashtab, attrib) in rdata
    v = env.value
    out = {}
    frame, hashtab = v.frame, v.hash_table
    if frame is not None and frame.info.type == T.LIST:
        out.update(pairlist(frame))
    if hashtab is not None and hashtab.info.type == T.VEC:
        for b in hashtab.value:
            if b is not None and b.info.type == T.LIST:
                out.update(pairlist(b))
    return out

def vec(o):
    o = deref(o)
    if o.info.type in (T.REAL, T.INT, T.LGL):
        a = np.ma.filled(np.ma.asarray(o.value).astype(float), np.nan)
        if o.info.type in (T.INT, T.LGL):
            a[a == -2147483648] = np.nan
        return a
    if o.info.type == T.STR:
        return strings(o)
    raise TypeError(o.info.type)

def dataframe(o):
    o = deref(o)
    at = attrs(o)
    cols = strings(at["names"])
    rn = deref(at["row.names"])
    rows = strings(rn) if rn.info.type == T.STR else None
    data = {}
    for c, v in zip(cols, o.value):
        v = deref(v)
        a = vec(v)
        va = attrs(v)
        if "levels" in va:  # factor
            lev = strings(va["levels"])
            a = [None if np.isnan(x) else lev[int(x) - 1] for x in a]
        data[c] = a
    return pd.DataFrame(data, index=rows)

def read_eset(path):
    p = rdata.parser.parse_file(path)
    top = pairlist(p.object)
    name, es = next(iter(top.items()))
    es = deref(es)
    slots = attrs(es)
    ad = env_vars(slots["assayData"])
    m = deref(ad["exprs"])
    ma = attrs(m)
    dim = [int(x) for x in m and ma["dim"].value]
    X = np.ma.filled(np.ma.asarray(m.value).astype(float), np.nan).reshape(dim[1], dim[0]).T
    dn = deref(ma["dimnames"]).value
    rows, cols = strings(dn[0]), strings(dn[1])
    expr = pd.DataFrame(X, index=rows, columns=cols)
    pheno = dataframe(attrs(deref(slots["phenoData"]))["data"])
    try:
        feat = dataframe(attrs(deref(slots["featureData"]))["data"])
    except Exception:
        feat = None
    return name, expr, pheno, feat
