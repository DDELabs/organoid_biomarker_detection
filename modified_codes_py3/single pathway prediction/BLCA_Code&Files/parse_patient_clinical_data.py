### parse patient clinical data ###
import scipy.stats as stat
import numpy as np
from collections import defaultdict
import pandas as pd
execfile('pathway_utilities.py', globals())
execfile('parse_Drugbank.py', globals())
gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()
common_dbID, dbID_commonSyn, commonSyn_dbID, common_Syn, Syn_common = parse_Drugbank_drugbankID_synonyms()


## ---------------------------------------------------------------
## patient data
def parse_TCGA_drug_treatment_data(cancer_type):
	'''
	output = { drug : [ list of patients ] }
	output2 = { patient : [drugs] }
	-----------------------------------------------
	drug synonyms are gathered from DRUGBANK database
	'''
	drugPat, patDrug = defaultdict(list), defaultdict(list) # { drug : [ list of patients ] }, { patient : [ drugs ] }
	if 'blca' in cancer_type.lower():
		fi_directory = './data/TCGA_BLCA/clinical_drug/9014e5dd-78fa-4507-bb07-156370362990/nationwidechildrens.org_clinical_drug_blca.txt'
	df = pd.read_table(fi_directory)
	for i in range(len(df)):
		pat, drug = df['bcr_patient_barcode'][i], df['pharmaceutical_therapy_drug_name'][i].upper().replace('-','').replace(' ','')
		if 'TCGA' in pat:
			drugs = drug.split('AND')
			for drug in drugs:
				if (not 'unknown'.upper() in drug.upper()) and (not 'notavailable'.upper() in drug.upper()):
					common_drug = drug
					if drug in Syn_common:
						common_drug = Syn_common[drug]
					drugPat[common_drug].append(pat)
					patDrug[pat].append(common_drug)
	return drugPat, patDrug

def parse_TCGA_survival_data_boolean_format(cancer_type):
	'''
	output = { patient : { 'survival_month' ; month, 'vital_status' : 1, 0 } }
	'''
	if 'blca' in cancer_type.lower():
		output = parse_TCGA_BLCA_survival_boolean_format()
	return output



## ---------------------------------------------------------------
## blca
def parse_TCGA_blca_drug():
	'''
	output = { drug : [ list of patients ] }
	output2 = { patient : [ drugs ] }
	--------------------------------------------------------
	drug names are converted to common drug names
	common drug names and drug synonyms are gathered from DRUGBANK database
	'''
	drugPat, patDrug = defaultdict(list), defaultdict(list) # { drug : [ list of patients ] }, { patient : [ drugs ] }
	fi_directory = './data/TCGA_BLCA/clinical_drug/9014e5dd-78fa-4507-bb07-156370362990/nationwidechildrens.org_clinical_drug_blca.txt'
	df = pd.read_table(fi_directory)
	for i in range(len(df)):
		pat, drug = df['bcr_patient_barcode'][i], df['pharmaceutical_therapy_drug_name'][i].upper().replace('-','').replace(' ','')
		if 'TCGA' in pat:
			common_drug = drug
			if drug in common_Syn:
				common_drug = drug
			else:
				if drug in Syn_common:
					common_drug = Syn_common[drug]
			drugPat[common_drug].append(pat)
			patDrug[pat].append(common_drug)
	return drugPat, patDrug



## patient survival
def parse_TCGA_BLCA_survival():
	'''
	output = { pat : { 'status' : 'alive'/'dead', 'months' : month } }
	'''
	output = {} # { pat : { 'status' : 'alive'/'dead', 'months' : month } }
	days_dic = defaultdict(list) # { 'days_to_death','days_to_last_followup' : [ days ] }

	fi_directory = './data/TCGA_BLCA/clinical_patient/eaa71705-960a-4abd-b5d7-f5fdc0d0c5af/nationwidechildrens.org_clinical_patient_blca.txt'
	df = pd.read_table(fi_directory)
	for i in range(len(df)):
		pat, status = df['bcr_patient_barcode'][i], df['vital_status'][i].upper()
		if ('TCGA' in pat):
			if status.upper() == 'ALIVE':
				days = df['last_contact_days_to'][i]
			elif status.upper() == 'DEAD':
				days = df['death_days_to'][i]

			# # check if days are correctly provided  
			if (not '[DISCREPANCY]' in days.upper()) and (not '[NOT AVAILABLE]' in days.upper()) and (not '[NOT APPLICABLE]' in days.upper()) and (not '[' in days):
				months = float(days)/(365./12.)
				if not pat in output:
					output[pat] = {}
				output[pat]['status'] = status
				output[pat]['months'] = months
	return output


def parse_TCGA_BLCA_survival_boolean_format():
	'''
	output = { patient ID : { 'survival_month'/'vital_status' } } // 1:deceased, 0:alive
	'''
	output = {}
	surDic = parse_TCGA_BLCA_survival()
	for pat in surDic:
		if ('months' in surDic[pat]) and ('status' in surDic[pat]):
			if not pat in output:
				output[pat] = {}
			if 'months' in surDic[pat]:
				output[pat]['months'] = surDic[pat]['months']
			if 'status' in surDic[pat]:
				if 'DEAD' in surDic[pat]['status']:
					output[pat]['status'] = 1
				elif 'ALIVE' in surDic[pat]['status']:
					output[pat]['status'] = 0
	result_output = {}
	for pat in output:
		if ('status' in output[pat]) and ('months' in output[pat]):
			result_output[pat] = output[pat]
	return result_output

