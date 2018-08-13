import matplotlib.pyplot as plt 
import os, sys
from rosetta import *
import pdb

from numpy import *
init(extra_options="-mute all")


def calculate_energy_pose(pose):
    global scorefxn
    global SS
    global iterCounter
    aux = 0
    sec_struct = Dssp(pose)
    sec_struct.insert_ss_into_pose(pose);
    SS2 = pose.secstruct()
    
    for j in range(len(SS)):
        if SS[j] == "H" or SS[j] == "G" or SS[j] == "I":
            if SS2[j] == "H":
                aux += -10.0
            else:
                aux += 10.0
        if SS[j] == "E" or SS[j] == "B" or SS[j] == "b":
            if SS2[j] == "E":
                aux += -10.0
            else:
                aux += 10.0
        if (SS[j] == "C" or SS[j] == "T") and SS2[j] != "L":
            aux += 10.0
    score = scorefxn(pose)
    iterCounter += 1
    #print SS2,score,aux
    #sasa = (calc_total_sasa(pose, 1.4))
    return score+aux#+sasa

rmsdDict = {#'1AB1':[],
'1K43':[],
'1ACW':[],
'1WQC':[],
'2MR9':[],

#'2P5K':[],
'2P81':[],
#'1AB1':[],
#'1CRN':[],
'1ROP':[],
'2KDL':[],
"2M7T":[],
"1ZDD":[],
"1UTG":[]}
#'2MTW':[]}

rmsdDictAdaptativo = {#'1AB1':[],
'1K43':[],
'1ACW':[],
'1WQC':[],
'2MR9':[],

#'2P5K':[],
'2P81':[],
#'1AB1':[],
#'1CRN':[],
'1ROP':[],
'2KDL':[],
"2M7T":[],
"1ZDD":[],
"1UTG":[]}
#'2MTW':[]}


# rmsdDictAdaptativo = {'1K43':[],
# '1ACW':[],
# #['1K43',"RGKWTYNGITYEGR","TTTEEETTEEECCC"],
# '1WQC':[],
# '2MR9':[],
# '2MTW':[],
# #['2P5K',"NKGQRHIKIREIITSNEIETQDELVDMLKQDGYKVTQATVSRDIKELHLVKVPTNNGSYKYSL","CHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHCCCCCHHHHHHHHHHHCCEEEEETTTEEEEEC"],
# '2P81':[],
# '1AB1':[]
# #['1CRN',"TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN","CEECCCHHHHHHHHHHHHCCCCHHHHHHHHCCEECCCCCCTTTTTC"],
# #['1ROP',"MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARF","CCHHHHHHHHHHHHHHHHHHHHHHHHHHHCCHHHHHHHHHHHHHHHHHHHHHHHHC"],
# #['2KDL',"TTYKLILNLKQAKEEAIKELVDAGTAEKYIKLIANAKTVEGVWTLKDEIKTFTVTE","TTTTTCCCHHHHHHHHHHHHHHHCCCHHHHHHHHCCCCHHHHHHHHHHHHHCCCCC"],
# #["2M7T","GCPQGRGDWAPTSCSQDSDCLAGCVCGPNGFCG","CCTTTBTTBTCEECCCGGGCTTTTCEETTTEEC"],
# #["1ZDD","FNMQCQRRFYEALHDPNLNEEQRNAKIKSIRDDC","CCHHHHHHHHHHHHTTTTTHHHHHHHHHHHHHHC"],
# #["1UTG","GICPRFAHVIENLLLGTPSSYETSLKEFEPDDTMKDAGMQMKKVLDSLPQTTRENIMKLTEKIVKSPLCM","CCCHHHHHHHHHHHHCCHHHHHHHHHHCCCCHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHCGGGC"]
# ]


# rmsdLowest = {'1K43':[],
# '1ACW':[],
# #['1K43',"RGKWTYNGITYEGR","TTTEEETTEEECCC"],
# '1WQC':[],
# '2MR9':[],
# '2MTW':[],
# #['2P5K',"NKGQRHIKIREIITSNEIETQDELVDMLKQDGYKVTQATVSRDIKELHLVKVPTNNGSYKYSL","CHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHCCCCCHHHHHHHHHHHCCEEEEETTTEEEEEC"],
# '2P81':[],
# '1AB1':[]
# #['1CRN',"TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN","CEECCCHHHHHHHHHHHHCCCCHHHHHHHHCCEECCCCCCTTTTTC"],
# #['1ROP',"MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARF","CCHHHHHHHHHHHHHHHHHHHHHHHHHHHCCHHHHHHHHHHHHHHHHHHHHHHHHC"],
# #['2KDL',"TTYKLILNLKQAKEEAIKELVDAGTAEKYIKLIANAKTVEGVWTLKDEIKTFTVTE","TTTTTCCCHHHHHHHHHHHHHHHCCCHHHHHHHHCCCCHHHHHHHHHHHHHCCCCC"],
# #["2M7T","GCPQGRGDWAPTSCSQDSDCLAGCVCGPNGFCG","CCTTTBTTBTCEECCCGGGCTTTTCEETTTEEC"],
# #["1ZDD","FNMQCQRRFYEALHDPNLNEEQRNAKIKSIRDDC","CCHHHHHHHHHHHHTTTTTHHHHHHHHHHHHHHC"],
# #["1UTG","GICPRFAHVIENLLLGTPSSYETSLKEFEPDDTMKDAGMQMKKVLDSLPQTTRENIMKLTEKIVKSPLCM","CCCHHHHHHHHHHHHCCHHHHHHHHHHCCCCHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHCGGGC"]
# ]

f = open("DictResultPSOCAN.txt", "w")

path = "/home/mariel/Desktop/SANPSO/"
pathC = "/home/mariel/Desktop/PSO_CAN/"
pathRef = "/home/mariel/Desktop/TestSet_PDBFiles/"
#getFiles = os.listdir(path)
for protein in rmsdDict.keys():
	print protein	
	#if protein != "1ACW": continue
	for i in range(30):
		pdbFile = protein+"_"+str(i)
		pose1 = pose_from_pdb(pathC+protein+"_"+str(i)+"/gBest.pdb")
		pose2 = pose_from_pdb(pathRef+protein+".pdb")
		rmsdDict[protein].append([CA_rmsd(pose1, pose2,2,pose2.n_residue()-2),i])
	for i in range(30):
		pdbFolder = protein+"_"+str(i)
		getFiles = os.listdir(path+pdbFolder)
		#print getFiles
		getRMSD = []
		for file in getFiles:
		# 	if "luster_1_" in file:
		# 		#print file
		# #try:
		# 		#pdbFile = "Cluster_1_"+pdbFolder
		# 		pose1 = pose_from_pdb(path+pdbFolder+"/"+file)
		# 		pose2 = pose_from_pdb(pathRef+protein+".pdb")
		# 		rmsdDictAdaptativo[protein].append(CA_rmsd(pose1, pose2,2,pose2.n_residue()-2))
		# rmsdDictAdaptativo[protein].append(min(getRMSD))
			if "luster_" in file:
				#print file
				pose1 = pose_from_pdb(path+pdbFolder+"/"+file)
				pose2 = pose_from_pdb(pathRef+protein+".pdb")
				rmsdDictAdaptativo[protein].append([CA_rmsd(pose1, pose2,2,pose2.n_residue()-2),file,i])

	f.write(str(protein)+'\n'+str(rmsdDict[protein])+"\n\n")



for k, v in rmsdDictAdaptativo.iteritems():
	#f.write(str(k)+'\n'+str(v)+"\n")
	print k,'\t', v




'''
#print rmsdLowest	
#pdb.set_trace()

getRuns = {'1AB1':[],
'1ACW':[],
'1K43':[],
'1WQC':[],
'2MR9':[],
'2MTW':[],
#'2P5K':[],
'2P81':[]}
getRuns2 = {'1AB1':[],
'1ACW':[],
'1K43':[],
'1WQC':[],
'2MR9':[],
'2MTW':[],
#'2P5K':[],
'2P81':[]}
for protein in rmsdDict.keys():
	meanProt = mean(rmsdDict[protein])
	stdProt = std((rmsdDict[protein]))
	getRuns2[protein] = [meanProt,stdProt]
	
f = open("stats24.02.txt", "w")
f.write("Proteina\tClusteNumber\tlowestRMSDfromAbundantCluster\n\n")
for k, v in getRuns.iteritems():
	f.write(str(k)+'\t'+str(v[0])+'\t'+str(v[1])+"\n")
# f.write("Proteina\tClusteNumber\tlowestRMSDofALL\n\n")
# for k, v in getRuns2.iteritems():
# 	f.write(str(k)+'\t'+str(v)+"\n")



'''
'''
except:
	pdbFile = "Bestfromcluster_1_"+pdbFolder
	pose1 = pose_from_pdb(path+pdbFolder+"/"+pdbFile+".pdb")
	pose2 = pose_from_pdb(pathRef+protein+".pdb")
	rmsdDictAdaptativo[protein].append(CA_rmsd(pose1, pose2,2,pose2.n_residue()-2))
		
		
		
		
xAxis = ['1AB1','1AB1_Nich','1ACW','1ACW_Nich','1K43','1K43_Nich','1WQC','1WQC_Nich', '2MR9','2MR9_Nich','2MTW','2MTW_Nich', '2P81','2P81_Nich' ]

fig = plt.figure(1, figsize=(12, 9))
ax = fig.add_subplot(111)
bp = ax.boxplot([rmsdDict['1AB1'],rmsdDictAdaptativo['1AB1'],rmsdDict['1ACW'],rmsdDictAdaptativo['1ACW'],rmsdDict['1K43'],rmsdDictAdaptativo['1K43'],rmsdDict['1WQC'],rmsdDictAdaptativo['1WQC'],rmsdDict['2MR9'],rmsdDictAdaptativo['2MR9'],rmsdDict['2MTW'],rmsdDictAdaptativo['2MTW'],rmsdDict['2P81'],rmsdDictAdaptativo['2P81']])
ax.set_xticklabels(xAxis,fontsize = 7)
fig.savefig('compFinaleiraMenorRMSD.png', bbox_inches='tight')

		'''
