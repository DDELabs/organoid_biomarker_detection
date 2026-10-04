# Fixing the cloud environment (network access + setup script)

## What went wrong

The domain list was pasted into the environment's **Setup script** field. That field is a
bash script run at the start of every new session, so each domain line runs as a command:

```
/tmp/init-script.sh: line 3: ftp.ncbi.nlm.nih.gov: command not found   (exit code 127)
```

As a result:
1. The domains were **never allowed**, because the network allowlist is a different setting.
2. **Every new session in this environment fails to start** until the setup script is fixed.
   (The current session was created before the edit, so it keeps running.)

## Fix, in about two minutes

Open the cloud environment menu in the session's title bar, then **Edit** on the "Default"
environment. Docs: https://code.claude.com/docs/en/cloud-environments#network-access

### 1. Setup script: replace everything with

```bash
#!/bin/bash
pip install -q pandas numpy scipy scikit-learn statsmodels gseapy networkx matplotlib joblib pytest rdata openpyxl zstandard || true
```

(Or leave it empty. It must contain only shell commands, never bare domain names.)

### 2. Network access: choose **Custom**, keep the default package-manager list, and add these under *Allowed domains*

```
api.gdc.cancer.gov
portal.gdc.cancer.gov
gdc.cancer.gov
docs.gdc.cancer.gov
ftp.ncbi.nlm.nih.gov
www.ncbi.nlm.nih.gov
eutils.ncbi.nlm.nih.gov
ncbi.nlm.nih.gov
stringdb-downloads.org
string-db.org
www.ebi.ac.uk
ftp.ebi.ac.uk
ega-archive.org
www.biosino.org
bioconductor.org
bioarchive.galaxyproject.org
depot.galaxyproject.org
zenodo.org
figshare.com
ndownloader.figshare.com
depmap.org
www.cancerrxgene.org
cog.sanger.ac.uk
cellmodelpassports.sanger.ac.uk
xenabrowser.net
synapse.org
www.synapse.org
repo-prod.prod.sagebase.org
www.cell.com
ars.els-cdn.com
www.nature.com
static-content.springer.com
europepmc.org
www.hgcc.se
ctrdb.ncpsb.org.cn
reactome.org
www.gsea-msigdb.org
data.broadinstitute.org
dgidb.org
huggingface.co
```

Alternatively, pick a broader access level (for example "Full") if your organisation allows it.

### 3. Start a new session

Network and setup-script changes apply only to sessions created **after** the edit. Ask the
new session to verify access:

```bash
for u in https://api.gdc.cancer.gov/status https://ftp.ncbi.nlm.nih.gov/geo/ https://stringdb-downloads.org/download/; do
  echo "$(curl -sS -o /dev/null -w '%{http_code}' --max-time 15 $u) $u"; done
```

Expect 200 (or 301/302/403-from-the-site) rather than `CONNECT tunnel failed, response 403`.

## Once access works

`robust/data_builders/README.md` lists what to fetch directly instead of from mirrors.
Most importantly:
- GSE72970 and the full GSE109211 (STORM, with its placebo arm);
- Lee 2018 bladder organoid drug data (Cell supplement plus GSE103990);
- GDC STAR counts via the API.
