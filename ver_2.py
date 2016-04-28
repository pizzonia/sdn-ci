#!/usr/bin/python

#
# I suggest to run under linux run with 
# stdbuf -o0 python ver_2.py | tee out.txt
#

from gurobipy import *



K = 1024
KBps = K*8 #1 kilobyte (per second) expressed in bit (per second)
MBps = K*KBps
GBps = K*K*KBps
Mbps = K*K # 1 megabit (per second)
Gbps = K*K*K # 1 gigabit (per second)



IDS='IDS'
nodes_N = [] # list of the network nodes
nodes_M = [IDS] # list of the device nodes

capacity = [] # map form pair of nodes (an arc) to capacity
arcs = [] # list of undirected arcs (they will be duplicated in directed arcs afterwards)


streams = []  # list of bidirectional streams (they will be duplicated in directed streams afterwards)
            # bidirectional streams, 6-tuple: node a, node b, bandwidth from a to b, bandwidth form b to a, relevance from a to b, relevance from b to a



####################################
# substation creation utilities

def transformationSystem(pfx):
    controlDevicesSingle= [  # follows wikipedia substation description
               #### primary 
               'pv0',# voltage meeter
               'ps1',# circuit switch
               'pb1',# breaker
               'pc1',# current transformer 
               'ps2',# circuit switch
               'pc2',# current transformer
               'pb2',# breaker
               ##### transformer
               'tr0',
               #### secondary
               'sb1',# breaker
               'sc1',# current transformer
               'ss1',# circuit switch
               'sc2',# current transformer
               'sb2',# breaker
               'ss2',# circuit switch
               'sv0'# voltage meeter
               ]
    controlDevicesSingle = map(lambda x: pfx+x, controlDevicesSingle)
    
    controlDevicesRedounded=[]
    for dev in controlDevicesSingle:
        controlDevicesRedounded.append(dev+'x')
        controlDevicesRedounded.append(dev+'y')
        controlDevicesRedounded.append(dev+'z')
        
    return controlDevicesRedounded
        

def fullBipariteGraph(devs1, devs2, capacity):
    "connect all devices of two sets (iterables)"
    from itertools import product
    return { (d1,d2): capacity for d1,d2 in product(devs1,devs2) }




####################################

def set_network_for_four_big_substations():
    # the network for four substations with three transformers each, fully redounded network
    global nodes_M, nodes_N, arcs, capacity, streams
    
    #mbw=10*KBps  # bandwidth demand for measurements
    mbw = 100*KBps  # bandwidth demand for measurements
    swcap = 1*Gbps # inter switch link capacity
    devcap = 10*Mbps # switch to embedded devices link capacity
    srvcap = 1*Gbps # switch to server link capacity
    sstcap = 4*Mbps # inter substation link capacities
    
    links={}

    switchesMainRoom = ['sw0', 'sw1']
    nodes_N += switchesMainRoom
    links[('sw0', 'sw1')] = swcap 
    links.update( fullBipariteGraph(switchesMainRoom, [IDS], srvcap)  )  

    for ss in ['N', 'S', 'E', 'W']:  # four substations

        switchesRoom = [ss+'sw3', ss+'sw4']
        nodes_N += switchesRoom
        links[(ss+'sw3',ss+'sw4')] = swcap 
        
        devicesRoom = [ss+'scada', ss+'hmi', ss+'historian']
        nodes_M += devicesRoom  
        links.update( fullBipariteGraph(switchesRoom, devicesRoom, swcap)  )
        
        for ts in ['A', 'B', 'C']:
#        for ts in ['A']:
            ctrlDevs = transformationSystem(ss+ts)
            nodes_M += ctrlDevs
                
            switchesField = [ss+ts+'sw1', ss+ts+'sw2']
            nodes_N += switchesField
            links[(ss+ts+'sw1',ss+ts+'sw2')] = swcap
            
            links.update( fullBipariteGraph(switchesField, ctrlDevs, devcap)  )  
            links.update( fullBipariteGraph(switchesRoom, switchesField, swcap)  )  

        ##### TO BE PUT IN THE PAPER AS A TABLE IN THE EVALUATION SECTION!!!!!!!
        # bidirectional streams, a, b, bandwidth from a to b, bandwidth form b to a, relevance from a to b, relevance from b to a
        streams += [ (ss+'scada',dev, mbw/10, mbw, 1, 100) for dev in ctrlDevs if dev[3]=='v' ]  # voltage meters
        streams += [ (ss+'scada',dev, mbw/10, mbw, 1, 100) for dev in ctrlDevs if dev[3]=='c' ]  # current meters
        streams += [ (ss+'scada',dev, 1.5*KBps, 1.5*KBps, 100, 10 ) for dev in ctrlDevs if dev[3]=='s' ]  # circuit switches
        streams += [ (ss+'scada',dev, 1.5*KBps, 1.5*KBps, 10, 100) for dev in ctrlDevs if dev[3]=='b' ]  # circuit breaker
        streams += [ (ss+'scada',dev, 5*mbw/10, 5*mbw, 50, 50) for dev in ctrlDevs if dev[2:3]=='tr' ]  # transformer
        streams += [ (ss+'scada',ss+'hmi', 300*mbw, 30*mbw, 50, 50) ]  
        streams += [ (ss+'scada',ss+'historian', 300*mbw, 30*mbw, 50, 50) ]  

    links.update( fullBipariteGraph(switchesMainRoom, ['Nsw3','Nsw4'], swcap)  )  
    
    for a,b in [('N','E'),('E','S'),('S','W'),('W','N')]:  # four inter substation links
        links[(a+'sw3',b+'sw4')] = sstcap  


    arcs, capacity = multidict(links)
    print "network topology and streams finished: "+ str(datetime.now())
    
    

def set_network_for_one_big_substations():
    # the network for a substation with three transformers, fully redounded network
    global nodes_M, nodes_N, arcs, capacity, streams

    #mbw=10*KBps  # bandwidth demand for measurements
    mbw = 100*KBps  # bandwidth demand for measurements
    swcap = 100*Mbps # inter switch link capacity
    devcap = 10*Mbps # switch to embedded devices link capacity
    srvcap = 1*Gbps # switch to server link capacity
    #sstcap = 4Mbps # inter substation link capacities
    
    
    links={}

    switchesRoom = ['sw3', 'sw4']
    nodes_N += switchesRoom
    links[('sw3','sw4')] = swcap 
    
    links.update( fullBipariteGraph(switchesRoom, [IDS], swcap)  )  
    
    devicesRoom = ['scada', 'hmi', 'historian']
    nodes_M += devicesRoom  
    links.update( fullBipariteGraph(switchesRoom, devicesRoom, swcap)  )
    
    for ts in ['A','B','C']:
        ctrlDevs = transformationSystem(ts)
        nodes_M += ctrlDevs
            
        switchesField = [ts+'sw1', ts+'sw2']
        nodes_N += switchesField
        links[(ts+'sw1',ts+'sw2')] = swcap
        
        links.update( fullBipariteGraph(switchesField, ctrlDevs, devcap)  )  
        links.update( fullBipariteGraph(switchesRoom, switchesField, swcap)  )  

    # bidirectional streams, a, b, bandwidth from a to b, bandwidth form b to a, relevance from a to b, relevance from b to a
    streams += [ ('scada',dev, mbw/10, mbw, 1, 100) for dev in ctrlDevs if dev[2]=='v' ]  # voltage meters
    streams += [ ('scada',dev, mbw/10, mbw, 1, 100) for dev in ctrlDevs if dev[2]=='c' ]  # current meters
    streams += [ ('scada',dev, 1.5*KBps, 1.5*KBps, 100, 10 ) for dev in ctrlDevs if dev[2]=='s' ]  # circuit switches
    streams += [ ('scada',dev, 1.5*KBps, 1.5*KBps, 10, 100) for dev in ctrlDevs if dev[2]=='b' ]  # circuit breaker
    streams += [ ('scada',dev, 5*mbw/10, 5*mbw, 50, 50) for dev in ctrlDevs if dev[1:2]=='tr' ]  # transformer
    streams += [ ('scada','hmi', 100*mbw, 10*mbw, 50, 50) ]  
    streams += [ ('scada','historian', 100*mbw, 10*mbw, 50, 50) ]  

    arcs, capacity = multidict(links)



def set_network_for_one_small_substation():
    # the network for a small substation with only one transformation system, fully redounded network
    global nodes_M, nodes_N, arcs, capacity, streams

    #mbw=10*KBps  # bandwidth demand for measurements
    mbw = 100*KBps  # bandwidth demand for measurements
    swcap = 10*Mbps # inter switch link capacity
    devcap = 10*Mbps # switch to embedded devices link capacity
    srvcap = 10*Mbps
    
    links={}
    
    ctrlDevs = transformationSystem("")
    nodes_M += ctrlDevs
        
    switchesField = ['sw1', 'sw2']
    nodes_N += switchesField
    links[('sw1','sw2')] = swcap
    
    links.update( fullBipariteGraph(switchesField, ctrlDevs, devcap)  )  
    
    switchesRoom = ['sw3', 'sw4']
    nodes_N += switchesRoom
    links[('sw3','sw4')] = swcap 

    links.update( fullBipariteGraph(switchesRoom, switchesField, swcap)  )  
    links.update( fullBipariteGraph(switchesRoom, [IDS], swcap)  )  

    devicesRoom = ['scada', 'hmi', 'historian']
    nodes_M += devicesRoom
    
    links.update( fullBipariteGraph(switchesRoom, devicesRoom, swcap)  )

    # bidirectional streams, a, b, bandwidth from a to b, bandwidth form b to a, relevance from a to b, relevance from b to a
    streams += [ ('scada',dev, mbw/10, mbw, 1, 100) for dev in ctrlDevs if dev[1]=='v' ]  # voltage meters
    streams += [ ('scada',dev, mbw/10, mbw, 1, 100) for dev in ctrlDevs if dev[1]=='c' ]  # current meters
    streams += [ ('scada',dev, 1.5*KBps, 1.5*KBps, 100, 10 ) for dev in ctrlDevs if dev[1]=='s' ]  # circuit switches
    streams += [ ('scada',dev, 1.5*KBps, 1.5*KBps, 10, 100) for dev in ctrlDevs if dev[1]=='b' ]  # circuit breaker
    streams += [ ('scada',dev, 5*mbw/10, 5*mbw, 50, 50) for dev in ctrlDevs if dev[0:1]=='tr' ]  # transformer

    arcs, capacity = multidict(links)


####################################
def set_medium_network():
    global nodes_M, nodes_N, arcs, capacity, streams

    #Model data
    nodes_N = ['sw1', 'sw2', 'sw3', 'sw4']  # list of switches
    nodes_M += ['n1', 'n2', 'n3', 'sc1'] # list of machines
    
    # indirect links
    arcs, capacity = multidict({
      ('sc1', 'sw1'):   100,
      ('sc1', 'sw2'):   100,
      ('sw1', 'sw2'):   100,
      ('sw1', 'sw3'):   100,
      ('sw1', 'sw4'):   100,
      ('sw2', 'sw3'):   100,
      ('sw2', 'sw4'):   100,
      ('sw3', 'sw4'):   100,
      ('n1', 'sw3'):   100,
      ('n1', 'sw4'):   100,
      ('n2', 'sw3'):   100,
      ('n2', 'sw4'):   100,
      ('n3', 'sw3'):   100,
      ('n3', 'sw4'):   100,
      ('IDS', 'sw1'):   100,
      ('IDS', 'sw2'):   100
      })
    
    
     # bidirectional streams, a, b, bandwidth from a to b, bandwidth form b to a, relevance from a to b, relevance from b to a
    streams = [
        ('n1', 'sc1', 50, 10, 1, 1 ),
        ('n2', 'sc1', 50, 10, 1, 1 ),
        ('n3', 'sc1', 50, 10, 1, 1 )
    ]



#######################################################
#  network specification

def set_small_network():
    global nodes_M, nodes_N, arcs, capacity, streams, relevance 
    nodes_N = ['sw1', 'sw2']  # list of switch
    nodes_M += ['n1', 'n2', 'n3', 'n4'] # list of machines
     
    arcs, capacity = multidict({
      ('sw1', 'sw2'):   100,
      ('n1', 'sw2'):   100,
      ('n2', 'sw2'):   100,
      ('n3', 'sw1'):   100,
      ('n4', 'sw1'):   100,
      ('sw1', 'IDS'):  100
    })
     
    # bidirectional streams, a, b, bandwidth from a to b, bandwidth form b to a, relevance from a to b, relevance from b to a
    streams = [
        ('n1', 'n3', 80, 10, 1, 1),
    ]



########################################

from datetime import datetime
print str(datetime.now()) + " start"

set_network_for_four_big_substations()
#set_network_for_one_big_substations()
#set_network_for_one_small_substation()
#set_medium_network()
#set_small_network()
    
    
# nodes specified in arcs should be only in nodes list
assert  set( reduce(lambda x,y:x+y, map(list,arcs))).issubset(set(nodes_N+nodes_M))

# nodes specified in streams should be only in nodes list
assert  set( reduce(lambda x,y:x+y, map(lambda x: [x[0],x[1]], streams) )).issubset(set(nodes_N+nodes_M))

assert IDS in set(nodes_M)


dirarcs = []  # directed arcs
c = {}  # capacity of directred links

for i,j in arcs:
    c [i,j] = capacity[i,j]
    c [j, i] = capacity[i,j]
    dirarcs.append( (i,j) )
    dirarcs.append( (j,i) )
dirarcs = tuplelist(dirarcs)


relevance = {}
dirstreams = [] # unidirectional streams (from, to, bandwidth, name)
for n1, n2, b_fw, b_rev, rel_fw, rel_rev in streams:
    sigmaFw = n1+'->'+n2
    sigmaRev = n2+'->'+n1
    dirstreams.append((n1, n2, b_fw, sigmaFw))
    dirstreams.append((n2, n1, b_rev, sigmaRev))
    relevance[sigmaFw]= rel_fw
    relevance[sigmaRev]= rel_rev


print str(datetime.now()) + " network created"
########################################
# helper

def neighors(v):
    return set([u for u, w in dirarcs.select('*', v) ]) | set([w for u, w in dirarcs.select(v, '*') ])




########################################
# GUROBI model setup 
m = Model('netflow')
m.setAttr("ModelSense", -1)   # maximization
m.params.timeLimit = 40.0 
m.params.NumericFocus = 3 
  
  
# GUROBI Variables creation

print str(datetime.now()) + " starting defining model"
print " vars: x and r..."
# x[sigma,s,t] flow of sigma for each directed edge s->t  
x = {}
r = {}
for _s,_t,b,sigma in dirstreams:
    for ss,tt in dirarcs:
        x[sigma,ss,tt] = m.addVar ( vtype=GRB.BINARY, obj=-float(b)/c[ss, tt], name='x_%s_%s_%s' % (sigma, ss, tt) )
        r[sigma,ss,tt] = m.addVar ( vtype=GRB.BINARY, obj=-float(b)/c[ss, tt], name='r_%s_%s_%s' % (sigma, ss, tt) )


print "expressions OutX, OutR, InX, InR, F, Fr"

OutX = {} #  how much stream actually exit from a vertex (gurobi expression)
for s,t,b,sigma in dirstreams:
    for v in nodes_N + nodes_M:
        OutX[sigma,v] = quicksum(x[ sigma,ss,tt] for ss, tt in dirarcs.select(v, '*'))

InX = {} #  how much stream actually enter in a vertex (gurobi expression)
for s,t,b,sigma in dirstreams:
    for v in nodes_N + nodes_M:
        InX[sigma,v] = quicksum(x[ sigma,ss,tt] for ss, tt in dirarcs.select('*',v))

OutR = {} #  how much replica stream actually exit from a vertex (gurobi expression)
for s,t,b,sigma in dirstreams:
    for v in nodes_N + nodes_M:
        OutR[sigma,v] = quicksum(r[ sigma,ss,tt] for ss, tt in dirarcs.select(v, '*'))

InR = {} #  how much replica stream actually enter in a vertex (gurobi expression)
for s,t,b,sigma in dirstreams:
    for v in nodes_N + nodes_M:
        InR[sigma,v] = quicksum(r[ sigma,ss,tt] for ss, tt in dirarcs.select('*',v))


F = {} # F[sigma, v] stream flow imbalance for each vertex and commodity (gurobi expression)
Fr = {} # Fr    replica stream flow imbalance
for _s,_t,_b,sigma in dirstreams:
    for v in nodes_M + nodes_N:
        F[sigma, v] = OutX[sigma,v] - InX[sigma,v] 
        Fr[sigma, v] = OutR[sigma,v] - InR[sigma,v] 

    
M=sum( 2*float(b)/c[ss,tt] for _s,_t,b,_sigma in dirstreams for ss,tt in dirarcs  )
print "M=", M

repl={}
for s,t,b,sigma in dirstreams:
    repl[sigma] = m.addVar ( vtype=GRB.BINARY, obj=M*relevance[sigma], name='isReplicated_%s' % (sigma) )

m.update()

print str(datetime.now()) + " variables defined"


# this is for correctly stating the objective function
for s,t,b,sigma in dirstreams:
    m.addConstr( repl[sigma] == InR[sigma, IDS] , name='replicated_%s' % (sigma) )



# constraints 

print str(datetime.now()) + " constr: capacity..."
# Arc capacity constraints
for ss, tt in dirarcs:
    m.addConstr(quicksum( b*(x[sigma,ss,tt]+r[sigma,ss,tt])  for _s,_t,b,sigma in dirstreams) 
                   <= c[ss,tt],
                'capacity_%s_%s' % (ss, tt))

print str(datetime.now()) + " constr: all critical streams should be routed..."
# impose all streams routed
for s,t,b,sigma in dirstreams:
    m.addConstr( OutX[sigma,s] == 1, 'routed_source_%s' % (sigma))
    m.addConstr( InX[sigma,t] == 1, 'routed_target_%s' % (sigma))


print str(datetime.now()) + " constr: flow conservation (normal)...."
# normal flow conservation
for s,t,_b,sigma in dirstreams:
    for v in set(nodes_N + nodes_M) - set([s,t]):
        m.addConstr( F[sigma, v] == 0, 'consX_imb_%s_%s' % (v, sigma))
    # normal flow source and sink
    #m.addConstr( F[sigma, s] == 1, 'sourceX_%s_%s' % (s, sigma))   # all flow routed, mandatory
    #m.addConstr( F[sigma, t] == -1, 'targetX_%s_%s' % (t, sigma))



# print "constr: In/Out flow <= 1..."
# # exiting and entering in a vertex: flow bounds for each stream
# for s,t,_b,sigma in dirstreams:
#     for v in set(nodes_N + nodes_M):
#         m.addConstr( OutX[sigma, v] <= 1, 'limitX_out_%s_%s' % (v, sigma))
#         m.addConstr( InX[sigma, v] <= 1, 'limitX_in_%s_%s' % (v, sigma))
#         m.addConstr( OutR[sigma, v] <= 1, 'limitR_out_%s_%s' % (v, sigma))
#         m.addConstr( InR[sigma, v] <= 1, 'limitR_in_%s_%s' % (v, sigma))


    # replica flow conservation
print str(datetime.now()) + " constr: replica flow conservation..."
for s,t,_b,sigma in dirstreams:
    lasthops = set( ss for ss, tt in dirarcs.select('*',t) )
    for v in set(nodes_N) - lasthops:
        m.addConstr( Fr[sigma, v] == 0, 'consR_%s_%s' % (v, sigma))    

    # replica stream starts from v only if v is the last hop of the selected path - see paper
print str(datetime.now()) + " constr: replicas starts only from last hop..."
for s,t,_b,sigma in dirstreams:
    lasthops = set( ss for ss, tt in dirarcs.select('*',t) )
    for v in lasthops: # for all possible last hop of sigma
        m.addConstr( x[sigma,v,t] >= Fr[sigma, v] 
                     , 'replicaIfLastHop_%s_%s' % (v, sigma))

print str(datetime.now()) + " constr: IDS is not a source..."
for s,t,_b,sigma in dirstreams:
    for ss, tt in dirarcs.select(IDS, '*'):
        m.addConstr( r[sigma,ss,tt]  == 0, 'replicaNoIDSloops_%s_%s' % (sigma,tt))




#     # replica: source equals target imbalance for replica streams
# print str(datetime.now()) + " constr: replica starts and end flow equals..."
# for s,t,_b,sigma in dirstreams:
#     lasthops = set( ss for ss, tt in dirarcs.select('*',t) )
#     m.addConstr( quicksum( Fr[sigma, v] for v in lasthops ) == - Fr[sigma,IDS], 
#                  'SrcEqTrgR_%s' % (sigma))

    # avoid exiting loops of replica streams from IDS





print str(datetime.now()) + " constr: M nodes cannot switch (normal)..."
for s,t,_b,sigma in dirstreams:
    # normal stream, nodes_M cannot switch
    for v in set(nodes_M) - set([s,t]):
        for ss,tt in dirarcs.select('*', v) + dirarcs.select(v, '*'):
            m.addConstr( x[sigma,ss,tt] == 0, 'MNodesNoSwitchX_%s_%s_%s_%s' % (v, sigma,ss,tt))

    # replica stream, nodes_M cannot switch
print str(datetime.now()) + " constr: M nodes cannot switch (replica)..."
for s,t,_b,sigma in dirstreams:
    for v in set(nodes_M) - set([IDS]):
        for ss,tt in dirarcs.select('*', v) + dirarcs.select(v, '*'):
            m.addConstr( r[sigma,ss,tt] == 0, 'MNodesNoSwitchR_%s_%s_%s_%s' % (v, sigma,ss,tt))

    





m.update()

print str(datetime.now()) + " model setup finished"

import os
print os.getcwd()
#m.write(os.path.join(os.getcwd(),'out.ilp') )
m.write('out.lp' )
# Compute optimal flow


print str(datetime.now()) + " starting optimization"

m.optimize()

print str(datetime.now()) + " optimization finished"


def recreateSequence(darcs):
    from toposort import toposort_flatten
    
    dependFrom={}
    for u,v in darcs:
        dependFrom[v]=set([u])

    return toposort_flatten(dependFrom)

# Print flow
if m.status in [ GRB.Status.OPTIMAL, GRB.Status.TIME_LIMIT ]  :
    flow = m.getAttr('x', x)
    replica = m.getAttr('x', r)
    isReplicated =  m.getAttr('x', repl)
    for s,t,_b,sigma in dirstreams:
        print('\nPath for stream %s:' % sigma)
        darc=[]
        for i,j in dirarcs:
            if flow[sigma,i,j] > 0:
                print('x %s -> %s: %g' % (i, j, flow[sigma,i,j]))
                darc.append((i,j))
        path =  recreateSequence(darc)
        assert path[0] == s
        assert path[-1] == t
        print path

        if not isReplicated[sigma]:
            continue
        
        darc=[]
        for i,j in dirarcs:
            if replica[sigma,i,j] > 0:
                print('r %s -> %s: %g' % (i, j, replica[sigma,i,j]))
                darc.append((i,j))
        replicapath = recreateSequence(darc)
        assert replicapath[-1] == IDS
        assert replicapath[0] == path[-2]
        print replicapath



    def cap(e):
        return c[e[0],e[1]]
    def occX(e):
        u=e[0]
        v=e[1]
        return sum( flow[sigma,u,v]*b for _s,_t,b,sigma in dirstreams)
    def occR(e):
        u=e[0]
        v=e[1]
        return sum( replica[sigma,u,v]*b for _s,_t,b,sigma in dirstreams)
    def occ(e):
        return occR(e)+occX(e)
    def occFrac(e):
        return float(occ(e))/cap(e)

    darcsByOccupation = [ (occFrac(e) ,e) for e in dirarcs ]
    darcsByOccupation.sort()
    for occ ,e in darcsByOccupation:
        if occ == 0: 
            continue
        print "%s->%s  cap=%d   occ=%d (%%%.3f)   (norm=%d , repl=%d)" % (e[0], e[1], cap(e), occX(e)+occR(e), float(occX(e)+occR(e))/cap(e)*100,  occX(e), occR(e) )

#m.printAttr("x")
#m.printAttr("NumConstrs")

