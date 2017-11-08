#model constructor
import copy
import igraph
import random
import copy
import pylab as pl
import scipy
from scipy import random
import heapq

class simpleNetworkSIRModel():
    def __init__(self, b = .2,  g = .01, S = 300, I = 1, p = .02, nei = 4):
        
        #parameters
        self.b = b  # The fixed number of contacts per day of each infected individual 
        self.g = g # 
        self.t = 0
        self.p = p

        self.N = S + I
        
        #make a small world graph with as many nodes as we have individuals
        self.graph = igraph.Graph.Watts_Strogatz(1, self.N, nei=nei, p = p)
        #we do this to get rid of multiple edges and self-loops that the 
        #randomly generated small-world graph might have
        self.graph.simplify()
        self.initgraph=copy.deepcopy(self.graph)
        
        
        #we're going to keep track of who is next to who using a list of lists
        self.adjacencyList = []
        for i in range(self.N):
            self.adjacencyList.append([])
        
        #now we're going to unpack the info from the graph 
        #into a more usable format
        
        #this is an efficient way of doing it bust shows you 
        #how to turn the graph into an 
        #adjacency list
        for edge in self.graph.es: #looping over the graph's edge sequence
            #indexing adjacency by node ID,so we can do quick lookups 
            self.adjacencyList[edge.source].append(edge.target) 
            self.adjacencyList[edge.target].append(edge.source) 
        
        #to use iGraph's internal method to do this, comment 
        #out above and uncomment the following:
        #self.adjacencyList = graph.get_adjlist()
                                                                
        #going to use this to store the *indices* of agents in each state
        self.sAgentList = []
        self.iAgentList = []
        self.rAgentList = []

        #and here we're going to store the counts of how many agents are in each
        #state @ each time step
        self.sList = []
        self.iList = []
        self.rList = []
        self.newIList = []
        
        #and we'll use this to keep track of recovery times in the more
        #efficient implementation
        self.recoveryTimesHeap = []

        #make a list of agent indices (easy because they're labeled 0-N)
        allAgents = [ i for i in range(self.N)]
        
        #shuffle the list so there's no accidental correlation in agent actions
        random.shuffle(allAgents)

        #start with everyone susceptible
        self.sAgentList = copy.copy(allAgents)
        
        #now infect a few to infect at t = 0
        self.indexCases = []
        for i in range(I):
            indexCase = self.sAgentList[0]
            self.indexCases.append(indexCase)
            self.infectAgent(indexCase)
            self.iAgentList.append(indexCase)
            
# heap-based method for recovering agents using an arbitrary distribution of recovery times

    def infectAgent(self,agent):
        self.sAgentList.remove(agent)
        
        #uncomment for exponentially distributed recovery times
        recoveryTime = self.t + scipy.random.exponential(1/self.g)
        
        #comment out above and uncomment below to try different recovery
        #distributions
        
        #lognormal with mean 1/g
        #recoveryTime = self.t + scipy.random.lognormal(mean = scipy.log(1/g), sigma = scipy.log(20) )
        #if recoveryTime <= self.t:
        #    recoveryTime = self.t + 1
        
        #gamma distributed times with mean 1/g
        #shape = 10.0
        #recoveryTime = self.t + scipy.random.gamma(shape, scale = 1/(shape*self.g))
        
        #normally distributed with mean 1/g
        #recoveryTime = self.t + scipy.random.normal(1/self.g, 10)
        #if recoveryTime <= self.t:
        #    recoveryTime = self.t + 1
        
        #constant
        #recoveryTime = self.t + (1/g)
        
        #note that we're pushing a tuple onto the heap where the first element
        #is the recovery time and the second one is the agent's unique ID
        heapq.heappush(self.recoveryTimesHeap, (recoveryTime, agent))
        
        return 1
    
    def recoverAgents(self):
        #when we recover agents, it's similar to the previous
        #non-network implementation
        recoverList = []
        if len(self.recoveryTimesHeap) > 0:
            while self.recoveryTimesHeap[0][0] <= self.t:
                #we take advantage of python's built-in sequence sorting methods
                #which compare starting from the first element in a sequence,
                #so if these are all unique, we can sort arbitary sequences 
                #by their first element without a special comparison operator
                recoveryTuple = heapq.heappop(self.recoveryTimesHeap)
                recoverList.append(recoveryTuple[1])
                if len(self.recoveryTimesHeap) == 0:
                    break
    
        return recoverList

           
    #again, the guts of the model
    def run(self):
        #same as while I > 0
        ctr=0
        while len(self.iAgentList) > 0:
            print("Executing While,iter",ctr)
            ctr+=1
            tempIAgentList = []
            recoverList = []
            newI = 0
            #we only need to loop over the agents who are currently infectious
            for iAgent in self.iAgentList:
                #and then expose their network neighbors
                for agent in self.adjacencyList[iAgent]:
                    #given that the neighbor is susceptible
                    if agent in self.sAgentList:
                        if (random.random() < self.b): 
                            #and then it's the same as the other models
                            newI += self.infectAgent(agent)
                            tempIAgentList.append(agent)
                             
            
            #then get the list of who is recovering
            recoverList = self.recoverAgents()
            
            #and do the bookkeeping with agent indices
            
            #for recoveries
            for recoverAgent in recoverList:
                self.iAgentList.remove(recoverAgent)
                self.rAgentList.append(recoverAgent)
            
            #and new infections
            self.iAgentList.extend(tempIAgentList)
            
            #then track the number of individuals in each state
            self.sList.append(len(self.sAgentList))
            self.iList.append(len(self.iAgentList))
            self.rList.append(len(self.sAgentList))
            self.newIList.append(newI)
            
            #increment the time
            self.t += 1
            #print('t', self.t, 'numS', len(self.sAgentList), 'numI', len(self.iAgentList) )
            
            #reshuffle the agent list so they step in a random order the next time
            #around
            random.shuffle(self.iAgentList)
            yield
        
        #and when we're done, return all of the relevant information
        return [self.sList, self.iList, self.rList, self.newIList]

    def graphPlotInit(self):
        for v in self.initgraph.vs():
            v['label_size']=0
        
            v['color'] = 'blue'
            if v.index in self.indexCases: #Initial infected nodes are green.
                v['color'] = 'green'
        
        if self.p <= .05:
            l = self.initgraph.layout_circle()
        
        #as p grows larger, it's probably more useful to use a spring-embedder
        #to see structure and clustering
        elif len(self.initgraph.vs) < 500: 
            #this is the most common spring-embdedder layout for graphs 
            #but can be inefficient for large graphs
            l = self.initgraph.layout_kamada_kawai()

        else:
            #but this one is faster and gives qualitatively
            #similar results
            l = self.initgraph.layout_grid_fruchterman_reingold()
        
        #now just plot the graph, passing it the layout we chose
        igraph.drawing.plot(self.initgraph)
        
    #this method lets us plot the network and the states of individuals at the end
    #of the run
    def graphPlot(self):
        for v in self.graph.vs():#loop over the nodes in the graph- in igraph the "vertex sequence"
            v['label_size'] = 0 #if you don't do this, nodes will be labeled with their indices, which
                                #can be visually confusing
                                
            v['color'] = 'blue' #all nodes start blue
#             if v.index in self.rAgentList or v.index in self.iAgentList:
#                 v['color'] = 'red' #if they were infected during the run, color them red
            if v.index in self.rAgentList:
                v['color'] = 'gray' #If they were infected and recovered.
            
            if v.index in self.iAgentList:
                v['color'] = 'red'
                
            if v.index in self.indexCases: #color index cases green
                v['color'] = 'green'
                
        # a circular layout can be the best when the disorder parameter is pretty low
        
        if self.p <= .05:
            l = self.graph.layout_circle()
        
        #as p grows larger, it's probably more useful to use a spring-embedder
        #to see structure and clustering
        elif len(self.graph.vs) < 500: 
            #this is the most common spring-embdedder layout for graphs 
            #but can be inefficient for large graphs
            l = self.graph.layout_kamada_kawai()

        else:
            #but this one is faster and gives qualitatively
            #similar results
            l = self.graph.layout_grid_fruchterman_reingold()
        
        #now just plot the graph, passing it the layout we chose
        igraph.drawing.plot(self.graph)
        #return self.graph