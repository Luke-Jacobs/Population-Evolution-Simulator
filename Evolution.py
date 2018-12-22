import random, matplotlib.pyplot as plt, numpy as np


class Trait:
    """Simulates a sequence of DNA that expresses itself in an Individual."""

    def __init__(self):
        """Default Trait properties."""
        self.name = ""
        self.genotype = []
        self.phenotype = None
        self.reproductiveFitness = None

    def __add__(self, other):
        """Combines Traits using shuffling during meiosis."""

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
        """For pretty printing."""
        return "<%s %s>" % (self.name, str(self.genotype))


class Individual:
    """A collection of traits, fitness, age, and a species identifier. Can mate with another Individual
    object to produce another Individual with reshuffled DNA (traits)."""

    def __init__(self):
        """Default Individual characteristics"""
        self.traits = []  # traits with different names
        self.speciesName = ""
        self.age = 0.0
        self.fitness = None

    def getTraitGenotype(self, name):
        """Lookup Trait object in Individual by name"""
        for trait in self.traits:
            if trait.name == name:
                return trait.genotype
        return None

    def __add__(self, other):
        """Mating with another individual."""

        offspring = Individual()

        if self.speciesName != other.speciesName:
            raise RuntimeError("Different species mating")

        #merge traits
        offspring.speciesName = self.speciesName
        for i in range(len(self.traits)):
            offspring.traits.append(self.traits[i] + other.traits[i])

        return offspring

    def getFitness(self, traitFitnesses):
        """Returns value of overall Individual fitness given a dictionary of [traitName: fitness]."""

        if self.fitness:  # If we already calculated this Individual's fitness
            return self.fitness

        self.fitness = 0.0
        for trait in self.traits:
            if trait.name in traitFitnesses: #if species that has trait possibility
                if trait.phenotype:
                    self.fitness += traitFitnesses[trait.name]
        return self.fitness

    def __str__(self):
        """For pretty printing."""
        traitMsg = ""
        for trait in self.traits:
            traitMsg += (str(trait) + " ")
        return "<%s %s>" % (self.speciesName, traitMsg)

    #make from genome (ordered list of trait names)
    #and from traitChances (list of percent change dominant allele)
    @staticmethod
    def fromGenome(speciesName: str, genome: list, traitChances=None):
        """Make an Individual from a genome and initial trait chances for an environment.
        @:param list genome: A list of trait names. ex: ['Red Eyes', 'Tall']
        @:param list traitChances: A list of dominant trait chances. ex: [.5, .9, 1.0]
        @:return Individual: The offspring of random chance."""

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
    """A collection of individuals."""

    def __init__(self, speciesInformation=None):
        self.individuals = []

        # Automatically populate this object
        if speciesInformation:
            self.addSpecies(*speciesInformation)

    def immigrate(self, individuals) -> None:
        self.individuals += individuals

    def show(self) -> None:
        """Prints the string representation of each individual in this population."""
        populationMsg = ""
        for individual in self.individuals:
            populationMsg += (str(individual) + "\n")
        print(populationMsg[:-1])

    def addSpecies(self, speciesName: str, traits: list, domAlleleChances: list, n: int) -> None:
        individuals = [Individual.fromGenome(speciesName, traits, domAlleleChances) for _ in range(n)]
        self.individuals += individuals

    def getAlleleFrequencies(self, traitNames) -> dict:
        """Return a dictionary of [traitName (str): frequency (float)] given the object's current state.
        Used in time-based simulation functions."""

        frequencyByTrait = {}  # [TraitName: frequency]
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

    def graphTraitFrequencies(self, traitNames: list) -> None:
        """Graph the number of traits in this population at one moment in time."""

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
            plt.ylabel('Number')
            plt.xticks(index, label)
            plt.title(traitName)
        plt.tight_layout()
        plt.show()


class Environment:
    """A simulated Environment that simplifies the effects of evolution on a population."""

    def __init__(self, population=None, traitDeathChances=None, traitMatingFitnesses=None):
        self.population = population if population else []
        self.traitDeathChances = traitDeathChances if traitDeathChances else {}  # Chances of insta-death for a trait
        self.traitMatingFitnesses = traitMatingFitnesses if traitMatingFitnesses else {}  # Benefits/disadvantages

    @staticmethod
    def traitDeath(population: Population, traitDeathChances: dict) -> Population:
        """Eliminates part of a population depending on each individual's traits."""

        dead = []
        # Mark individuals for death
        for individual in population.individuals:
            for trait in individual.traits:
                if trait.name in traitDeathChances:
                    if trait.phenotype:  # If genes are expressed (if one allele is dominant)
                        deathChance = traitDeathChances[trait.name]
                        if random.random() < deathChance:  # Individual dies
                            dead.append(individual)
        # Kill each death-marked individual
        for deadindividual in dead:
            try:
                population.individuals.remove(deadindividual)
            except ValueError:
                continue
        # Return the smaller population
        return population

    def mateTheMostFit(self, n=1, freqMate=.5) -> None:
        """Perform selection using constants that dictate random chance."""

        offspring = [] #new offspring
        individualsRanked = [] #individuals to mate

        #rank individuals based on traits
        for i in range(len(self.population.individuals)):
            individualsRanked.append([i, self.population.individuals[i].getFitness(self.traitMatingFitnesses)])
        individualsRanked.sort(key=lambda x: -x[1])

        # mate a part of the best individuals
        nWhoWillMate = int(len(self.population.individuals) * freqMate)
        # put individuals into list
        bestIndividuals = [self.population.individuals[ranking[0]] for ranking in individualsRanked]
        for i in range(0, nWhoWillMate, 2):
            try:
                pair = [bestIndividuals[i], bestIndividuals[i+1]]
            except KeyError:  # if odd list
                continue
            for j in range(n):  # have n babies
                offspring.append(pair[0] + pair[1]) #one offspring

        self.population.individuals += offspring

    def fastforward(self, timesteps: int):
        """Move forward a certain number of timesteps."""
        for _ in range(timesteps):
            self.mateTheMostFit(3, .5)
            self.population = self.traitDeath(self.population, self.traitDeathChances)


class Analysis:
    """A class for storing analytical functions in a convenient namespace."""

    @staticmethod
    def plotAlleleFrequencies(traitData: list):
        """Plots a line graph of allele frequencies over time given a list of trait data.
        :param list traitData: A list of dictionaries [traitname (str): traitfrequency (float)]"""

        graphPoints = [[] for _ in range(len(traitData[0]))] #should be list of separate trait points
        key = [traitName for traitName in traitData[0]]

        for timepoint in traitData: #for each time period
            for i in range(len(timepoint)): #for each trait number
                traitName = list(timepoint.keys())[i]
                graphPoints[i].append(timepoint[traitName])

        x = np.arange(len(traitData))
        for traitPoints in graphPoints:
            plt.plot(x, traitPoints)
        plt.legend(key)
        plt.xlabel('Generations')
        plt.ylabel('Frequency of allele')
        plt.title("Allele frequency vs. Time")
        plt.show()


if __name__ == "__main__":
    frogTraits = ["Crazy Color", "Long Tongue", "Green Eyes"]
    domChances = [.9, .5, .1]

    pop = Population(("Frog", frogTraits, domChances, 1000))
    pop.graphTraitFrequencies(frogTraits)

    env = Environment(population=pop)
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
    Analysis.plotAlleleFrequencies(data)

