import numpy as np
import mdtraj as md

ind=np.zeros((5001,40))
ind2=np.zeros((5001,40))
xxx=np.zeros((5001,40))
nwat=np.zeros(5001)
file1=open('index_shared_molecules.dat','r')  # to do the test
###file1=open('index_shared_molecules.dat','r')
k=-1
kll=-1
lt=0
for line in file1:
    lt=lt+1
    #print(lt)
    if(line[0:5]=="MODEL"):
        if(k>-1):
            nwat[k]=kll+1
        k=k+1
        kll=-1
    else:
        kll=kll+1
        #print(line)
        ind[k,kll]=int(line[0:7])
        xxx[k,kll]=float(line[6:30])


nwat[k]=kll+1
#print(nwat)
#print(xxx)
#print(ind)

file1.close()


file2=open('traj_wat.pdb','r')
km=-1
for line in file2:
    #print(line)
    if(line[0:5]=="MODEL"):
        km=km+1
        kcheck=0
        ij=0
        ind2[km,ij]=int(nwat[km])
        indd=0
        #print('here',km)
    if(line[0:4]=="ATOM"):
        #print(line)
        indd=indd+1
        name=line[13:15]
        if(ij < int(nwat[km])):
            if(name=='OW'):
    
                res=int(line[20:26])
                index=int(line[4:13])
                xmd=float(line[27:38])
                #print('hh',index)
                for ll in range(0,int(nwat[km])):
                    if(res==int(ind[km,ll]) and xmd==xxx[km,ll]):
                        #print('here')
                        #print(line)
                        #print(line)
                        #print(index)
                        #print(indd)
                        #print(int(ind[km,ll]))
                        ij=ij+1
                        ind2[km,ij]=int(indd)-1 #start from 0 mdtraj
                        #print(ij,nwat[km])
                        #stop
            #if(ij==int(nwat[km])):
                #print('out')
            #    break
      
print(nwat[0])

#index wat per frame ind2
print('wat',ind2[0:2,0:(int(nwat[0])+1)])


name=input("Enter name trajctory: ");
topname=input("Enter name topology with extension: ");
traj = md.load(name+'.xtc',top=topname)

frames=traj.n_frames
print(frames)

top=traj.topology
id_prot=top.select("protein and (element O or element P or element C or element N or element S)")

print('index protein')
print(id_prot)

CUTOFF=0.4

filewat=open('list_pairs_wat_prot.dat','w')
for time in range(0,frames):
    heavy_pairs = np.array(
        [(i,j) for  i in id_prot  for j in ind2[time,1:(int(nwat[time])+1)]])
    print(heavy_pairs)
    stringa='MODEL '+str(time)+'\n'
    filewat.write(stringa)
    heavy_pairs_distances_traj = md.compute_distances(traj[time], heavy_pairs)
    print('here')
    new_ind=heavy_pairs[heavy_pairs_distances_traj[0,:] < CUTOFF]
    print(new_ind)
    nel=len(new_ind)
    for i in range(0,nel):
        a1='%8d'%(int(new_ind[i][0])+1)
        a2='%8d'%(int(new_ind[i][1])+1)
        aresp='%8d'%(top.atom(int(new_ind[i][0])).residue.index+1)
        aresw='%8d'%(top.atom(int(new_ind[i][1])).residue.index+1)
        stringa=a1+'  '+a2+'  '+aresp+'  '+aresw+'\n'
        filewat.write(stringa)
        
        
        #       print(new_ind[i][0])
        #       stop
    #print(heavy_pairs_distances_traj)
    #print(len(heavy_pairs_distances_traj[0,:]),len(heavy_pairs))
    #print(heavy_pairs_distances_traj[0,:])
    #print(np.unique(heavy_pairs[heavy_pairs_distances_traj[0,:] < CUTOFF][:,0]))
    #new_ind=(np.unique(heavy_pairs[heavy_pairs_distances_traj[0,:] < CUTOFF][:,0]))
    #print(heavy_pairs[int(new_ind[0])])
    #print(heavy_pairs_distances_traj[0,int(new_ind[0])])
    #print(heavy_pairs[666])
    #print(heavy_pairs[0])
    #print(heavy_pairs[int(new_ind[0])])

    heavy_pairs=[]
    heavy_pairs_distances_traj=[]
    new_ind=[]

