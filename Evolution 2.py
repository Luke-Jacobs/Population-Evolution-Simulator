import random, matplotlib.pyplot as plt, numpy as np

class Trait:
    def __init__(self):
        self.name = ""
        self.genotype = []
        self.phenotype = None
        self.reproductiveFitness = None
    def __add__(self, other):
        if other.name != self.name:
            raise RuntimeError("Different traits mixing")

        newTrait = Trait()
        newTrait.name = self.name
        al1 = random.choice(self.genotype) #1 allele from self
        al2 = random.choice(other.genotype) #1 allele from other
        newTrait.genotype = [al1, al2]
        newTrait.phenotype = max(newTrait.genotype)
        return newTrait
    def __str__(self):
        return "<%s %s>" % (self.name, str(self.genotype))

class Individual:
    def __init__(self):
        self.traits = [] #traits with different names
        self.speciesName = ""
        self.age = 0.0
        self.fitness = None

    def getTraitGenotype(self, name):
        for trait in self.traits:
            if trait.name == name:
                return trait.genotype
        return None

    #mating with another individual
    def __add__(self, other):
        offspring = Individual()

        if self.speciesName != other.speciesName:
            raise RuntimeError("Different species mating")

        #merge traits
        offspring.speciesName = self.speciesName
        for i in range(len(self.traits)):
            offspring.traits.append(self.traits[i] + other.traits[i])

        return offspring

    #sets value of self.fitness given a dictionary of traitName: fitness
    def getFitness(self, traitFitnesses):
        if self.fitness:
            return self.fitness

        self.fitness = 0.0
        for trait in self.traits:
            if trait.name in traitFitnesses: #if species that has trait possibility
                if trait.phenotype:
                    self.fitness += traitFitnesses[trait.name]
        return self.fitness

    #for printing
    def __str__(self):
        traitMsg = ""
        for trait in self.traits:
            traitMsg += (str(trait) + " ")
        return "<%s %s>" % (self.speciesName, traitMsg)

    #make from genome (ordered list of trait names)
    #and from traitChances (list of percent change dominant allele)
    @staticmethod
    def fromGenome(speciesName: str, genome: list, traitChances=None):
        offspring = Individual()
        offspring.speciesName = speciesName

        if not traitChances:
            traitChances = [.5]*len(genome) #halvsies
        for i in range(len(genome)):
            newTrait = Trait()
            newTrait.name = genome[i]
            al1 = random.random() < traitChances[i]
            al2 = random.random() < traitChances[i]
            newTrait.genotype = [al1, al2]
            newTrait.phenotype = max(newTrait.genotype)
            offspring.traits.append(newTrait)

        return offspring

class Population:
    def __init__(self):
        self.individuals = []

    def immigrate(self, individuals) -> None:
        self.individuals += individuals

    # def mateAll(self, n=1) -> None:
    #     for i in range(n):
    #         #new offspring
    #         offspring = []
    #
    #         #random mating
    #         random.shuffle(self.individuals)
    #         pairs = [self.individuals[i:i+2] for i in range(0, len(self.individuals), 2)]
    #             #put individuals into pairs
    #
    #         #combine genes and add individuals to a new array
    #         for pair in pairs:
    #             if len(pair) != 2:
    #                 continue
    #             offspring.append(pair[0] + pair[1])
    #
    #         self.individuals += offspring

    def show(self):
        populationMsg = ""
        for individual in self.individuals:
            populationMsg += (str(individual) + "\n")
        print(populationMsg[:-1])

    def addSpecies(self, speciesName, traits, domAlleleChances, n):
        individuals = [Individual.fromGenome(speciesName, traits, domAlleleChances) for i in range(n)]
        self.individuals += individuals

    def getAlleleFrequencies(self, traitNames):
        frequencyByTrait = {} #traitName: frequency
        totalAlleles = len(self.individuals) * 2

        for i in range(len(traitNames)):
            traitName = traitNames[i] #specific trait to plot
            homoDom = homoRec = hetero = 0

            for individual in self.individuals:
                genotype = individual.getTraitGenotype(traitName)
                if genotype[0] == genotype[1]: #homozygous
                    if genotype[0]:
                        homoDom += 1
                    else:
                        homoRec += 1
                else: #heterozygous
                    hetero += 1
            frequencyByTrait[traitName] = (homoDom*2 + hetero) / totalAlleles
        return frequencyByTrait

    def graphTraitFrequencies(self, traitNames):
        plt.figure(1)
        for i in range(len(traitNames)):
            traitName = traitNames[i] #specific trait to plot
            plt.subplot(len(traitNames), 1, i + 1)
            homoDom = homoRec = hetero = 0

            for individual in self.individuals:
                genotype = individual.getTraitGenotype(traitName)
                if genotype[0] == genotype[1]: #homozygous
                    if genotype[0]:
                        homoDom += 1
                    else:
                        homoRec += 1
                else: #heterozygous
                    hetero += 1

            label = ["Homozygous Dominant", "Heterozygous", "Homozygous Recessive"]
            index = np.arange(len(label))
            plt.bar(index, [homoDom, hetero, homoRec])
            plt.ylabel('Frequency')
            plt.xticks(index, label)
            plt.title(traitName)
        plt.tight_layout()
        plt.show()

class Environment:
    def __init__(self):
        self.population = None
        self.traitDeathChances = {} #chances of insta-death for a specific trait
        self.traitMatingFitnesses = {} #benefits/disadvantages to reproductive fitness

    @staticmethod
    #eliminates part of a population depending on each individual's traits
    #traitDeathChances = {"blue color": .9, "tongue size": .1}
    def traitDeath(population: Population, traitDeathChances: dict):
        dead = []
        for individual in population.individuals:
            for trait in individual.traits:
                if trait.name in traitDeathChances:
                    if trait.phenotype: #if genes are expressed
                        deathChance = traitDeathChances[trait.name]
                        if random.random() < deathChance: #individual dies
                            dead.append(individual)
        for deadindividual in dead:
            try:
                population.individuals.remove(deadindividual)
            except ValueError:
                continue
        return population

    def mateTheMostFit(self, n=1, freqMate=.5) -> None:
        offspring = [] #new offspring
        individualsRanked = [] #individuals to mate

        #rank individuals based on traits
        for i in range(len(self.population.individuals)):
            individualsRanked.append([i, self.population.individuals[i].getFitness(self.traitMatingFitnesses)])
        individualsRanked.sort(key=lambda x: -x[1])

        #mate a part of the best individuals
        nWhoWillMate = int(len(self.population.individuals) * freqMate)
        #put individuals into list
        bestIndividuals = [self.population.individuals[ranking[0]] for ranking in individualsRanked]
        for i in range(0, nWhoWillMate, 2):
            try:
                pair = [bestIndividuals[i], bestIndividuals[i+1]]
            except KeyError: #if odd list
                continue
            for j in range(n): #have n babies
                offspring.append(pair[0] + pair[1]) #one offspring

        self.population.individuals += offspring

    #move forward x timesteps
    def fastforward(self, timesteps: int):
        for i in range(timesteps):
            self.mateTheMostFit(3, .5)
            self.population = self.traitDeath(self.population, self.traitDeathChances)

class Analysis:

    @staticmethod
    def plotAlleleFrequencies(traitData: list):
        graphPoints = [[] for i in range(len(traitData[0]))] #should be list of separate trait points
        key = [traitName for traitName in traitData[0]]

        for timepoint in traitData: #for each time period
            for i in range(len(timepoint)): #for each trait number
                traitName = list(timepoint.keys())[i]
                graphPoints[i].append(timepoint[traitName])

        x = np.arange(len(traitData))
        for traitPoints in graphPoints:
            plt.plot(x, traitPoints)
        plt.legend(key)
        plt.title("Allele frequency vs. Time")
        plt.show()

if __name__ == "__main__":
    frogTraits = ["Crazy Color", "Long Tongue", "Green Eyes"]
    domChances = [.9, .5, .1]

    pop = Population()
    pop.addSpecies("Frog", frogTraits, domChances, 1000)
    pop.graphTraitFrequencies(frogTraits)

    env = Environment()
    env.population = pop
    env.traitDeathChances = {"Crazy Color": .0,
                             "Long Tongue": .0,
                             "Green Eyes": .0}
    env.traitMatingFitnesses = {"Crazy Color": 0,
                             "Long Tongue": 0,
                             "Green Eyes": 0}

    generations = 10
    data = []
    for i in range(generations):
        env.fastforward(1)
        timestop = env.population.getAlleleFrequencies(frogTraits)
        data.append(timestop)
        print(timestop)
    Analysis.plotAlleleFrequencies(data)

