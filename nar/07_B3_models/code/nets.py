"""Neural models (torch, CPU): MOLI-style deep classifier (paper Fig 3f comparator), drug-fingerprint MLP
(pretrained on cell lines, fine-tuned on lab models), and a graph convolutional network (GCN) over STRING.
"""
import copy
from functools import lru_cache

import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn

from b3common import paper_network, targets

torch.set_num_threads(2)


def _seed(s):
    torch.manual_seed(s)
    np.random.seed(s % (2 ** 31))


# ================================================================================ MOLI-style deep network (Fig 3f)
class Moli(nn.Module):
    """Sharifi-Noghabi et al. 2019 style: encoder -> embedding, classifier on embedding; loss = BCE + gamma x triplet."""

    def __init__(self, d, h, drop):
        super().__init__()
        self.enc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(drop))
        self.cls = nn.Sequential(nn.Dropout(drop), nn.Linear(h, 1))

    def forward(self, x):
        e = self.enc(x)
        return self.cls(e).squeeze(-1), nn.functional.normalize(e, dim=1)


def _triplet(emb, y, margin):
    pos, neg = torch.where(y == 1)[0], torch.where(y == 0)[0]
    if len(pos) < 2 or len(neg) < 1:
        return emb.sum() * 0
    d = torch.cdist(emb, emb) ** 2
    same = (y[:, None] == y[None, :]) & ~torch.eye(len(y), dtype=bool)
    diff = y[:, None] != y[None, :]
    hardest_pos = torch.where(same, d, torch.full_like(d, -1)).max(1).values
    hardest_neg = torch.where(diff, d, torch.full_like(d, 1e9)).min(1).values
    ok = (hardest_pos >= 0) & (hardest_neg < 1e9)
    return torch.relu(hardest_pos - hardest_neg + margin)[ok].mean() if ok.any() else emb.sum() * 0


def moli_fit(X, yb, hp, seed=0):
    _seed(seed)
    X = torch.tensor(X, dtype=torch.float32)
    y = torch.tensor(yb, dtype=torch.float32)
    m = Moli(X.shape[1], hp["nodes"], hp["dropout"])
    opt = torch.optim.Adagrad(m.parameters(), lr=hp["lr"], weight_decay=hp["weight_decay"])
    n, bs = len(y), hp["batch_size"]
    for ep in range(hp["epochs"]):
        m.train()
        perm = torch.randperm(n)
        for i in range(0, n, bs):
            b = perm[i:i + bs]
            if len(b) < 2:
                continue
            logit, emb = m(X[b])
            loss = nn.functional.binary_cross_entropy_with_logits(logit, y[b]) + hp["gamma"] * _triplet(emb, y[b], hp["margin"])
            opt.zero_grad()
            loss.backward()
            opt.step()
    m.eval()
    return m


def moli_predict(m, X):
    with torch.no_grad():
        return torch.sigmoid(m(torch.tensor(X, dtype=torch.float32))[0]).numpy()


MOLI_GRID = [dict(batch_size=8, nodes=n, lr=lr, dropout=dr, weight_decay=1e-3, gamma=g, epochs=60, margin=mg)
             for n in (32, 128) for lr in (0.01, 0.05) for dr in (0.3, 0.6) for g in (0.1, 0.5) for mg in (0.5, 1.5)]


# ================================================================================ drug-fingerprint MLP
class FPNet(nn.Module):
    def __init__(self, p, h=128, hd=32):
        super().__init__()
        self.cell = nn.Sequential(nn.Linear(p, h), nn.ReLU(), nn.Dropout(0.2), nn.Linear(h, 64), nn.ReLU())
        self.drug = nn.Sequential(nn.Linear(p, hd), nn.ReLU())
        self.head = nn.Sequential(nn.Linear(64 + hd, 64), nn.ReLU(), nn.Linear(64, 1))

    def forward(self, xc, xd):
        return self.head(torch.cat([self.cell(xc), self.drug(xd)], 1)).squeeze(-1)


def fp_train(Xc, Fd, y, epochs=10, lr=1e-3, seed=0, bs=256, net=None, wd=1e-5):
    _seed(seed)
    Xc, Fd, y = (torch.tensor(a, dtype=torch.float32) for a in (Xc, Fd, y))
    net = net or FPNet(Xc.shape[1])
    opt = torch.optim.Adam(net.parameters(), lr=lr, weight_decay=wd)
    n = len(y)
    for ep in range(epochs):
        net.train()
        perm = torch.randperm(n)
        for i in range(0, n, bs):
            b = perm[i:i + bs]
            loss = ((net(Xc[b], Fd[b]) - y[b]) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
    net.eval()
    return net


def fp_predict(net, Xc, fd):
    with torch.no_grad():
        Xc = torch.tensor(Xc, dtype=torch.float32)
        Fd = torch.tensor(np.broadcast_to(fd, Xc.shape) if np.ndim(fd) == 1 else fd, dtype=torch.float32)
        return net(Xc, Fd).numpy()


def fp_finetune(net0, Xc, fd, y, steps=100, lr=1e-3, lam=1e-2, seed=0):
    """Fine-tune all weights on lab models, L2 pull towards the pretrained weights (keeps cell-line knowledge)."""
    _seed(seed)
    net = copy.deepcopy(net0)
    ref = [p.detach().clone() for p in net.parameters()]
    Xc = torch.tensor(Xc, dtype=torch.float32)
    Fd = torch.tensor(np.broadcast_to(fd, Xc.shape).copy(), dtype=torch.float32)
    y = torch.tensor(y, dtype=torch.float32)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    net.train()
    for _ in range(steps):
        loss = ((net(Xc, Fd) - y) ** 2).mean() + lam * sum(((p - r) ** 2).sum() for p, r in zip(net.parameters(), ref))
        opt.zero_grad()
        loss.backward()
        opt.step()
    net.eval()
    return net


# ================================================================================ GCN over STRING
@lru_cache(None)
def _net_arrays():
    g = paper_network()
    nodes = list(g.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    rows, cols = zip(*[(idx[a], idx[b]) for a, b in g.edges])
    A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(nodes),) * 2)
    A = (A + A.T).tocsr()
    deg = np.asarray(A.sum(0)).ravel()
    W = A @ sp.diags(1 / np.maximum(deg, 1))          # column-normalised transition matrix
    return nodes, idx, A, W


@lru_cache(None)
def rwr(drug, restart=0.5, iters=50):
    """Random walk with restart from the drug's targets on the paper's STRING>700 network: node -> visit probability."""
    nodes, idx, A, W = _net_arrays()
    t = [idx[x] for x in targets(drug) if x in idx]
    p0 = np.zeros(len(nodes))
    if not t:
        return dict(zip(nodes, p0))
    p0[t] = 1 / len(t)
    p = p0.copy()
    for _ in range(iters):
        p = (1 - restart) * (W @ p) + restart * p0
    return dict(zip(nodes, p))


def gcn_graph(drug, genes, var_order, n_prop=300, n_var=1000):
    """Node set: top-n_prop RWR genes + top-n_var most variable genes, restricted to available genes; induced subgraph."""
    nodes, idx, A, W = _net_arrays()
    avail = set(genes) & set(idx)
    pr = rwr(drug)
    top = [g for g, _ in sorted(pr.items(), key=lambda kv: -kv[1]) if g in avail][:n_prop]
    var = [g for g in var_order if g in avail][:n_var]
    keep = list(dict.fromkeys(top + var))
    ii = [idx[g] for g in keep]
    S = A[ii][:, ii] + sp.eye(len(ii))
    d = np.asarray(S.sum(1)).ravel()
    S = sp.diags(d ** -0.5) @ S @ sp.diags(d ** -0.5)
    S = S.tocoo()
    At = torch.sparse_coo_tensor(np.vstack([S.row, S.col]), S.data, S.shape, dtype=torch.float32).coalesce()
    prior = np.log10(np.array([pr[g] for g in keep]) + 1e-8)
    prior = (prior - prior.mean()) / (prior.std() or 1)
    return keep, At, torch.tensor(prior, dtype=torch.float32)


class GCN(nn.Module):
    def __init__(self, A, prior, h=16):
        super().__init__()
        self.A, self.prior = A, prior
        self.l1, self.l2 = nn.Linear(2, h), nn.Linear(h, h)
        self.att = nn.Parameter(torch.zeros(1))
        self.out = nn.Linear(2 * h, 1)

    def prop(self, H):                      # H: (B, N, c) -> A @ H over nodes
        B, N, c = H.shape
        Z = torch.sparse.mm(self.A, H.permute(1, 0, 2).reshape(N, B * c))
        return Z.reshape(N, B, c).permute(1, 0, 2)

    def forward(self, X):                   # X: (B, N) standardised expression
        B, N = X.shape
        H = torch.stack([X, self.prior.expand(B, N)], -1)
        H = torch.relu(self.prop(self.l1(H)))
        H = torch.relu(self.prop(self.l2(H)))
        w = torch.softmax(self.prior * torch.exp(self.att), 0)
        r = torch.cat([H.mean(1), (H * w[None, :, None]).sum(1)], 1)
        return self.out(r).squeeze(-1)


def gcn_fit(A, prior, X, y, steps=250, lr=0.01, wd=1e-4, seed=0):
    _seed(seed)
    m = GCN(A, prior)
    X, y = torch.tensor(X, dtype=torch.float32), torch.tensor(y, dtype=torch.float32)
    opt = torch.optim.Adam(m.parameters(), lr=lr, weight_decay=wd)
    bs = 256
    for _ in range(steps):
        b = torch.randperm(len(y))[:bs] if len(y) > bs else slice(None)
        loss = ((m(X[b]) - y[b]) ** 2).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()
    m.eval()
    return m


def gcn_predict(m, X):
    with torch.no_grad():
        return m(torch.tensor(X, dtype=torch.float32)).numpy()
