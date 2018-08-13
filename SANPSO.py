
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

def useHistogram(prob_list):
    global chi_boo
    prob_radius = 0.5
    floatprecision = 3
    id_probs = 0
    rnd = random.random()
    i = 0
    for w in prob_list:
        rnd -= w[2]
        if rnd < 0:
            id_probs = i
            break
        i += 1
    probs = prob_list[id_probs]
    phi = probs[0]+random.uniform(-prob_radius, prob_radius)
    psi = probs[1]+random.uniform(-prob_radius, prob_radius)
    angles = random.choice(probs[3])
    omega = angles[0]+random.uniform(-prob_radius, prob_radius)
    if chi_boo:
        chi1 = 999.9
        chi2 = 999.9
        chi3 = 999.9
        chi4 = 999.9
        if angles[1] != 999.9:
            chi1 = angles[1]+random.uniform(-prob_radius, prob_radius)
        if angles[2] != 999.9:
            chi2 = angles[2]+random.uniform(-prob_radius, prob_radius)
        if angles[3] != 999.9:
            chi3 = angles[3]+random.uniform(-prob_radius, prob_radius)
        if angles[4] != 999.9:
            chi4 = angles[4]+random.uniform(-prob_radius, prob_radius)

        aa_angles = [round(phi,floatprecision), round(psi, floatprecision), round(omega, floatprecision), round(chi1, floatprecision), round(chi2, floatprecision), round(chi3, floatprecision), round(chi4, floatprecision)] #angles[5]
    else:
        aa_angles = [round(phi,floatprecision), round(psi, floatprecision), round(omega, floatprecision)]
    return aa_angles

def carregaAPL(heat_maps, chi_boo, protein, AA, SS):

  pathGlobalAPL = "/home/mariel/Dropbox/TCC/"+protein+"/"

  pathGlobal_saida ="/home/mariel/Dropbox/TCC/PSO/"
  
  list_apls = []
  if heat_maps == 3:
    list_apls = [protein+"-3", protein+"-2-Right", protein+"-2-Left", protein+"-1"]
  if heat_maps == 1:
    list_apls = [protein+"-1"]

  print ["{0: >2}".format(i) for i in AA]
  print ["{0: >2}".format(i) for i in SS]
  print ["{0: >2}".format(str(i)) for i in range(len(AA))]
  print list_apls

  secondary_structure_sequence = copy.deepcopy(SS)
  #print secondary_structure_sequence
  amino_acid_sequence = copy.deepcopy(AA)
  #print amino_acid_sequence

  dirs_apl = [x for x in os.walk(os.path.join(str(pathGlobalAPL)+str(),'.')).next()[1] if not x.startswith(".")]

  for folder in dirs_apl:
    if folder.startswith(protein[:4]):
        if folder not in APL.keys():
            APL[folder] = {}
            if folder not in list_apls: continue
            files_apl = os.listdir(pathGlobalAPL+folder)
            for file_apl in files_apl:
                id_apl = '_'.join(file_apl.split("_")[:-1])
                if file_apl[-4:] == ".dat" and file_apl[-4:] not in APL[folder].keys():
                    APL[folder][id_apl] = []
                    with open(pathGlobalAPL+folder+"/"+file_apl, 'r') as arq:
                        for line in arq:
                            list_line = line.split()
                            APL[folder][id_apl].append([float(list_line[0]), float(list_line[1]), float(list_line[2]), eval(''.join(list_line[3:]))])
                    APL[folder][id_apl].sort(key = lambda x: x[2], reverse=True)
                    if len(APL[folder][id_apl]) == 0:
                        sai_daqui = APL[folder].pop(id_apl)
  return APL

def makeAPL_Pop(heat_maps, tamanhoPop, APL, amino_acid_sequence, secondary_structure_sequence, protein):
  estruturas = []
  for i in range(tamanhoPop):
    estrutura = {}
    for j in range(len(amino_acid_sequence)):
        estrutura[j] = []
        aa_angles = []
        if j != 0 and j != len(amino_acid_sequence)-1 and heat_maps == 3:
            AA_Ant = amino_acid_sequence[j-1]
            AA = amino_acid_sequence[j]
            AA_Prox = amino_acid_sequence[j+1]
            SS_Ant = secondary_structure_sequence[j-1]
            SS = secondary_structure_sequence[j]
            SS_Prox = secondary_structure_sequence[j+1]
            chance = random.random()
            prob = []
            if 0.0 <= chance < 0.5:
                key = AA_Ant+AA+AA_Prox+"_"+SS_Ant+SS+SS_Prox
                if key in APL[protein+"-3"]:
                    prob = APL[protein+"-3"][key]
                else:
                    if int(random.random()):
                        key = AA+AA_Prox+"_"+SS+SS_Prox
                        if key in APL[protein+"-2-Right"]:
                            prob = APL[protein+"-2-Right"][key]
                        else:
                            key = AA+"_"+SS
                            prob = APL[protein+"-1"][key]
                    else:
                        key = AA_Ant+AA+"_"+SS_Ant+SS
                        if key in APL[protein+"-2-Left"]:
                            prob = APL[protein+"-2-Left"][key]
                        else:
                            key = AA+"_"+SS
                            prob = APL[protein+"-1"][key]
            elif 0.5 <= chance < 0.75:
                if int(random.random()):
                    key = AA+AA_Prox+"_"+SS+SS_Prox
                    if key in APL[protein+"-2-Right"]:
                        prob = APL[protein+"-2-Right"][key]
                    else:
                        key = AA+"_"+SS
                        prob = APL[protein+"-1"][key]
                else:
                    key = AA_Ant+AA+"_"+SS_Ant+SS
                    if key in APL[protein+"-2-Left"]:
                        prob = APL[protein+"-2-Left"][key]
                    else:
                        key = AA+"_"+SS
                        prob = APL[protein+"-1"][key]
            else:
                key = AA+"_"+SS
                prob = APL[protein+"-1"][key]
            aa_angles = useHistogram(prob)
            estrutura[j] = aa_angles
        else:
            prob = []
            AA = amino_acid_sequence[j]
            SS = secondary_structure_sequence[j]
            key = AA+"_"+SS
            prob = APL[protein+"-1"][key]
            aa_angles = useHistogram(prob)
            estrutura[j] = aa_angles
    estruturas.append(estrutura)
  return estruturas


from math import *
def SA(gbest, w_temp = 1):
    global AA
    global SS
    best_energy=gbest[1]
    global energy
    print "Melhor inicial SA: "+str(best_energy), gbest[1]
    pose = pose_from_sequence(AA, 'centroid')
    for j in range(len(AA)):
        pose.set_phi(j+1, gbest[0][j][0])
        pose.set_psi(j+1, gbest[0][j][1])
        #pose.set_omega(j+1, gbest[0][j][2])
        #for kk in range(1,pose.residue(j+1).nchi()+1):
         #   try:
          #      pose.set_chi(kk, j+1,gbest[0][j][2+kk])
           # except:
            #    break;    
   # print "teste: "+str(calculate_energy_pose(pose))           
    radius=0.5*w_temp
    amino_order = list(AA)
    energy_curr = best_energy
    for i in range(len(AA)):
        bit_start = -radius
        bit_end = radius
        temperature = 1000*w_temp    
        max_succ=1
        total_it=1
        best_phi = gbest[0][i][0]
        backup_phi = gbest[0][i][0]
        phi = gbest[0][i][0]
        while(temperature > 0.1):
            cont_int=0
            succ_int=0
            #loop interno p/ msma temp
            while (cont_int < total_it) and (succ_int<max_succ):
                bit = random.uniform(bit_start, bit_end)
                #phi_ant = phi #valor anterior de phi
                phi += bit
                if (phi > (backup_phi - radius)) and (phi < (backup_phi + radius)):
                    #alteracao temporaria
                    pose.set_phi(i+1, phi)
                    #self.pocket_list[0].calculate_energy()
                    energy_temp = calculate_energy_pose(pose)
                    #print "temp: "+str(energy_temp)
                    #print "global: "+str(best_energy)
                    if(energy_temp < energy_curr) or (exp((energy_curr-energy_temp)/temperature) > random.random()):
                        energy_curr = energy_temp
                        if energy_temp < best_energy:
                            best_energy = energy_temp
                            succ_int+=1
                            best_phi = phi
                            #print "1", phi
                            print "SA minimizou phi--> "+str(best_energy)
                            #pose.dump_pdb(pathGlobal+str(protein)+"_PSO_"+str(strike)+"/Minimizados/"+str(best_energy)+".pdb")
                    else:
                        pose.set_phi(i+1,gbest[0][i][0])
                else:
                    pose.set_phi(i+1,gbest[0][i][0])
        
                cont_int+=1
            temperature = temperature*0.85 #0.85 #0.98
        gbest[0][i][0] = best_phi
        pose.set_phi(i+1,best_phi)
        #print "2", gbest[0][i][0]
        bit_start = -radius
        bit_end = radius
        temperature = 1000*w_temp    
        max_succ=1
        total_it=1
        best_psi = gbest[0][i][1]
        backup_psi = gbest[0][i][1]
        psi = gbest[0][i][1]
        energy_curr = best_energy
       
        while(temperature > 0.1):
            cont_int=0
            succ_int=0
            #loop interno p/ msma temp
            while (cont_int < total_it) and (succ_int<max_succ):
                bit = random.uniform(bit_start, bit_end)
                #psi_ant = psi #valor anterior de phi
                psi += bit
                if (psi > (backup_psi - radius)) and (psi < (backup_psi + radius)):
                    #alteracao temporaria
                    pose.set_psi(i+1, psi)
                    #self.pocket_list[0].calculate_energy()
                    energy_temp = calculate_energy_pose(pose)
                    #print "temp: "+str(energy_temp)
                    #print "global: "+str(best_energy)
                    if(energy_temp < energy_curr) or (exp((energy_curr-energy_temp)/temperature) > random.random()):
                        energy_curr = energy_temp
                        if energy_temp < best_energy:
                            best_energy = energy_temp
                            succ_int+=1
                            best_psi = psi
                            print "SA minimizou psi--> "+str(best_energy)
                            #pose.dump_pdb(pathGlobal+str(protein)+"_PSO_"+str(strike)+"/Minimizados/"+str(best_energy)+".pdb")
                    else:
                        pose.set_psi(i+1,gbest[0][i][1])
                else:
                    pose.set_psi(i+1,gbest[0][i][1])
        
                cont_int+=1
            temperature = temperature*0.85 #0.85 #0.98
        gbest[0][i][1] = best_psi
        pose.set_psi(i+1,best_psi)

    gbest[1] = best_energy

    pose=Pose()
    pose = pose_from_sequence(AA, 'centroid')
    for j in range(len(AA)):
        pose.set_phi(j+1, gbest[0][j][0])
        pose.set_psi(j+1, gbest[0][j][1])

    #print calculate_energy_pose(pose)

    print "Melhor depois do SA: "+str(gbest[1]), best_energy
    return gbest        

def PSO(particles,gbest,ind):
    iter_max = 100
    #iter_max = 10
    c1 = 4
    c2 = 2
    w1 = 0.9
    w2 = 0.4
    w = 0
    #init_time = datetime.datetime.now()
    #gbest = particles[0]
    k = 0
    breaker = 0
    convergenceFlag = 0
    localCounter = 0
    print "Swarm size ", len(particles)
    while k<iter_max:
        if (localCounter >= 1000):
            break;
        counter = 0
        print "-->"+ str(k)
        for p in particles:

            amino_sequence = AA
            pose=Pose()
            pose = pose_from_sequence(amino_sequence, 'centroid')
            for j in range(len(amino_sequence)):
                pose.set_phi(j+1, p[0][j][0])
                pose.set_psi(j+1, p[0][j][1])
                pose.set_omega(j+1,p[0][j][2])

            fitness =  calculate_energy_pose(pose)
            localCounter += 1 
            totalFitness.append(fitness)
            w = w + (w2-w1)*(float(k)/iter_max)
            #c2 = c1_1*gni + c1_2
            if fitness < p[1]: #Individual Best
                p[1] = fitness
                p[3] = p[0]

            if fitness < gbest[1]:#Global Best
                gbest = p
                loops = 0
                print "minimizou: "+str(fitness)
                breaker = -1
                convergenceFlag = -10
            
                 
                #pose.dump_pdb(pathGlobal+str(protein)+"_"+str(strike)+"/Minimizados/"+"Cluster"+str(ind)+"_"+str(fitness)+".pdb")
            for i in range(len(p[0])):
                for j in range(len(p[0][i])):
                    if len(p[0][i]) != len(gbest[0][i]):
                        continue
                    c11 = p[3][i][j] - p[0][i][j]
                    c22 = gbest[0][i][j] - p[0][i][j]
                    v =  w*p[2] + (c1 * random.random() * c11) + (c2 * random.random() * c22)    
                    p[0][i][j] = p[0][i][j] + v
                    #print "Velocity: "+str(v)
            totalFitness[counter] = p[1]
            counter += 1
        k += 1
        breaker += 1  
        if breaker == 100:
            break;
    return particles,gbest,convergenceFlag


import numpy
import random
from rosetta import *
init(extra_options="-mute all")
import pdb
import os
import copy
from rosetta.core.scoring.dssp import Dssp
import datetime
import time
import scipy.spatial.distance as ssd


import numpy as np
#import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

iterCounter = 0

APL = {}

SS = sys.argv[3]
AA = sys.argv[2]
strike = sys.argv[4]
APL = {}
protein = sys.argv[1]
iterCounter = 0



chi_boo = 0 # sem chi
#chi_boo = 1 # com chi
heat_maps = 3 #apl 3 2 e 1
#heat_maps = 1 #apl 1
APL = carregaAPL(heat_maps, chi_boo,protein, AA, SS)

population = []


tamanhoPop = 1000
#tamanhoPop = 50

population = makeAPL_Pop(heat_maps, tamanhoPop, APL, AA, SS, protein)

scorefxn = create_score_function('score3')
pathGlobal = "/home/mariel/Desktop/TCC/SANPSO/" 
os.system("mkdir "+pathGlobal+""+str(protein)+"_"+str(strike)+"/;")
os.system("mkdir "+pathGlobal+""+str(protein)+"_"+str(strike)+"/Inicial/;")
os.system("mkdir "+pathGlobal+""+str(protein)+"_"+str(strike)+"/Final/;")

class Particle:
    pass

particles = []
for i in range(len(population)):
    p = {} 
    p[0] = population[i].values() #angles
    p[1] = 11111111111110.0 #fitness
    p[2] = 0.0 #velocity
    p[3] = p[0] #individual best
    p[4] = []
    particles.append(p)

def hClustering(particles,AA):
    print "Clustering..."
    flag = 0
    distanceMatrix = []
    beg = 1
    # let the first particle be the global best
    totalFitness = []
    for p in range(len(particles)-1):
        pose=Pose()
        pose = pose_from_sequence(AA, 'centroid')
        for j in range(len(AA)):
            pose.set_phi(j+1, particles[p][0][j][0])
            pose.set_psi(j+1, particles[p][0][j][1])
            pose.set_omega(j+1,particles[p][0][j][2])
        for p2 in range(beg,len(particles)):
            pose2=Pose()
            pose2 = pose_from_sequence(AA, 'centroid')
            for j in range(len(AA)):
                pose2.set_phi(j+1, particles[p2][0][j][0])
                pose2.set_psi(j+1, particles[p2][0][j][1])
                pose2.set_omega(j+1,particles[p2][0][j][2])
            #print "calculating", p, p2
            particles[p][4].append(CA_rmsd(pose, pose2,2,pose2.n_residue()-2))
            distanceMatrix.append(CA_rmsd(pose, pose2,2,pose2.n_residue()-2))
        #print "calculating", p
        beg += 1

    data_link =  linkage(distanceMatrix, method = "complete")
    threshold = len(AA)**(1.0/3)
    numNiches = len(AA)
    while True:
        if numNiches <= 100:
            break;

        clusters = fcluster(data_link,t = threshold,criterion = 'distance')
        numNiches = len(list(set(clusters)))
        print "-- n niches",numNiches
        print "++ threshold",threshold

        threshold += .1
    return clusters


clusters = hClustering(particles,AA)
uniqueClusters = list(set(clusters))

swarm = dict((uniqueClusters[k], []) for k in range(len(uniqueClusters)))
for i in range(len(particles)):
        swarm[clusters[i]].append(particles[i])


print "Initial number of clusters from APL:", len(swarm.keys())


totalFitness = []
ind = 0
for p in particles:
    amino_sequence = AA
    pose=Pose()
    pose = pose_from_sequence(amino_sequence, 'centroid')
    for j in range(len(amino_sequence)):
        pose.set_phi(j+1, p[0][j][0])
        pose.set_psi(j+1, p[0][j][1])
        pose.set_omega(j+1,p[0][j][2])
#        for kk in range(1,pose.residue(j+1).nchi()+1):
#            try:
#                pose.set_chi(kk, j+1,p[0][j][2+kk])
#            except:
#                break;           
    fitness =  calculate_energy_pose(pose)
    p[1] = fitness
    pose.dump_pdb(pathGlobal+str(protein)+"_"+str(strike)+"/Inicial/"+str(ind)+".pdb")
    ind += 1
    totalFitness.append(fitness)
menorInicial = min(totalFitness)
minEnergy = min(particles, key=lambda x:x[1])

gbest = particles[particles.index(minEnergy)]

print "Minimal energy on cluster", clusters[particles.index(minEnergy)], ":", menorInicial

print "************* Beginning NichingPSO"
iteration = 0
nFlags = 0
flag = -100
while (iterCounter!=1000000):
#while (iterCounter<20):
    ind = 0
    newSwarm = []
    if len(swarm.keys())> 1:
        for subswarm in swarm.keys():
            ind += 1
            print "======> SubSwarm ",ind
            print "++ Convergence Flag:",nFlags
            if len(swarm[subswarm]) < 2:
                print "One particle swarm"
                subSA = SA(swarm[subswarm])
                newSwarm += subSA
                continue 
            localMin = min(swarm[subswarm], key=lambda x:x[1])
            subS,lbest, flag = PSO(swarm[subswarm],localMin,ind)
             
            newSwarm += subS
            #break;
            clusters = hClustering(newSwarm,AA)
            uniqueClusters = list(set(clusters))
            swarm = dict((uniqueClusters[k], []) for k in range(len(uniqueClusters)))
            for i in range(len(particles)):
		        swarm[clusters[i]].append(particles[i])
    else:
        break;


c = 1        
for subswarm in swarm.keys():
    pose=Pose()
    amino_sequence = AA
    lbestEnergy = min(swarm[subswarm], key=lambda x:x[1])
    lbest = swarm[subswarm][swarm[subswarm].index(lbestEnergy)]
    print "SA for cluster",subswarm
    lbest = SA(lbest)
    pose = pose_from_sequence(amino_sequence, 'centroid')
    for j in range(len(AA)):
        pose.set_phi(j+1, lbest[0][j][0])
        pose.set_psi(j+1, lbest[0][j][1])
        pose.set_omega(j+1,lbest[0][j][2])       
    pose.dump_pdb(pathGlobal+""+str(protein)+"_"+str(strike)+"/Bestfromcluster_"+str(c)+"_"+str(protein)+"_"+str(strike)+".pdb")
    c += 1 
it = 0
for subswarm in swarm.keys():
	for p in swarm[subswarm]:
		pose=Pose()
		amino_sequence = AA
		pose = pose_from_sequence(amino_sequence, 'centroid')
		for j in range(len(amino_sequence)):
		    pose.set_phi(j+1, p[0][j][0])
		    pose.set_psi(j+1, p[0][j][1])
		    pose.set_omega(j+1,p[0][j][2])       
		pose.dump_pdb(pathGlobal+str(protein)+"_"+str(strike)+"/Final/"+str(it)+".pdb")
		it += 1 


f = open(""+pathGlobal+""+str(protein)+"_"+str(strike)+"/results.txt", "w")
f.write('\nParticle Swarm Optimisation\n')
f.write( str(protein)+'\n')
f.write( 'ACABOU\n'+'-'*9)

f.write( '\nRESULTS\n'+'-'*7)
#f.write( '\titerations: '+str(k)+'\n')
f.write('-'*9)
f.close()
'''

#data_dist = pdist(distanceMatrix)
data_link = linkage(distanceMatrix)

print data_link

test = fcluster(data_link,t = 2.0,criterion = 'distance')

dendrogram(test,labels=range(len(particles)))
plt.xlabel('Samples')
plt.ylabel('Distance')
plt.suptitle('Samples clustering', fontweight='bold', fontsize=14)
plt.xlabel('Samples')
plt.ylabel('Distance')
plt.suptitle('Samples clustering', fontweight='bold', fontsize=14);
plt.savefig("fClust.png")
# from fastcluster import *
# %timeit data_link = linkage(data_array, method='single', metric='euclidean', preserve_input=True)
# dendrogram(data_link,labels=range(len(particles)))
# plt.xlabel('Samples')
# plt.ylabel('Distance')
# plt.suptitle('Samples clustering', fontweight='bold', fontsize=14);
# plt.savefig("fastcluster.png"')
'''
