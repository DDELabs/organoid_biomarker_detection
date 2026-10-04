# Data security checklist for controlled-access data

These are typical requirements in EGA DAAs (UMCU, ICR, HKU, Hartwig) and in the NIH Genomic Data Sharing security best practices. Tick each item, fill in the specifics, and attach the list to requests. Items marked **(confirm)** vary by DAC, so check the specific DAA.

## Storage and encryption
- [ ] Data is stored only on institution-managed servers/HPC in [country/EU]. It is never stored on personal laptops, USB drives or personal cloud (Dropbox, Google Drive).
- [ ] Encryption at rest: full-disk or volume encryption (AES-256, e.g. LUKS) on the storage holding raw files.
- [ ] Encryption in transit: downloads only via pyega3, the EGA Live Outbox, or the SRA Toolkit (TLS). Internal transfers via SSH/SFTP only.
- [ ] Crypt4GH / dbGaP `.ngc` keys are kept with owner-only permissions (chmod 600). They are never committed to git.
- [ ] Data is kept outside any git repository; add `data/controlled/` to `.gitignore`. No raw or individual-level files go to GitHub, Hugging Face or Zenodo.

## Access control
- [ ] Named users only, matching the approved application. Each user has an individual account; no shared logins.
- [ ] Multi-factor authentication for remote access. VPN required off-site.
- [ ] Unix group / ACL restricts the dataset directory to approved users (`chmod 750`, group = project).
- [ ] Access is removed within [5] working days when a user leaves. The DAC is notified of personnel changes **(confirm)**.
- [ ] Access audit logging is enabled.

## Cloud and third parties
- [ ] No commercial cloud (AWS, GCP, Azure, Colab) without **prior written DAC approval** **(confirm)**. For dbGaP, cloud use must be declared in the DAR.
- [ ] If cloud is approved: an EU/approved region, encrypted buckets, no public links, and a provider DPA signed by the institution.
- [ ] Individual-level data is not uploaded to external web tools or LLM/AI services.

## System hygiene
- [ ] Patched OS, firewall and endpoint protection. Backups are encrypted and access-controlled. Derived individual-level files (count matrices, VCFs, per-patient tables) are protected like raw data.

## Governance and GDPR
- [ ] Ethics/IRB approval: [IRB/ethics number]. The institution's DPO has been informed. A GDPR legal basis is documented, and for non-EU recipients an international transfer clause is in place **(confirm)**.
- [ ] No re-identification attempts and no linkage to other individual-level sources.
- [ ] Only aggregate results are published. Small-cell suppression applies (for example, no per-patient tables with rare clinical attributes).
- [ ] Breach procedure: notify the institution's security office and the DAC within [72 h].
- [ ] Annual or progress reports are filed if required (dbGaP: annual renewal).
- [ ] At project end, all copies are securely deleted (`shred`/crypto-erase) and a **written destruction certificate** is sent to the DAC.

**Responsible contacts:** PI [name]; IT/security officer [name]; institutional signatory [name].
