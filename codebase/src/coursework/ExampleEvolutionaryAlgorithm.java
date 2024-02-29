package coursework;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import model.Fitness;
import model.Individual;
import model.NeuralNetwork;
import model.StringIO;
public class ExampleEvolutionaryAlgorithm extends NeuralNetwork {
	@Override
	public void run() {
		//Initialise a population of Individuals with 0 weights
		System.out.println(Parameters.printParams());
		population = initialiseAll0();
		//Record a copy of the best Individual in the population
		best = population.get(getBest()).copy();
		System.out.println("Best From Initialisation " + best);
		SimpleDateFormat format = new SimpleDateFormat("MM-dd");
		SimpleDateFormat formaat = new SimpleDateFormat("HH-mm");
		Date date = new Date();
		String filePrefix = format.format(date);
		String time = formaat.format(date);
		StringIO.writeStringToFile(filePrefix + ".csv", "\n" + time + ", " + String.format("%.3f", best.fitness) + ", ", true);
		while (evaluations < Parameters.maxEvaluations) {

			// Select 2 Individuals from the current population.
			int p1idx = selectTournament();
			int p2idx = selectTournament();
			Individual parent1 = population.get(p1idx).copy();
			while(p2idx == p1idx){ // no same p's
				p2idx = selectTournament();
			}
			Individual parent2 = population.get(p2idx).copy();

			// Generate a child by crossover
			ArrayList<Individual> children = crossoverUniform(parent1, parent2);

			//mutate the offspring
			mutate(children);

			// Evaluate the children
			evaluateIndividuals(children);

			// Replace children in population
			replaceIfBetter(children);

			// check to see if the best has improved
			best = population.get(getBest()).copy();

			// Implemented in NN class.
			outputStats();

			//Increment number of completed generations
		}

		StringIO.writeStringToFile(filePrefix + ".csv", String.format("%.3f", best.fitness) + ", ", true);

		//save the trained network to disk
		saveNeuralNetwork();
	}


	 // Sets the fitness of the individuals passed as parameters (whole population)
	private void evaluateIndividuals(ArrayList<Individual> individuals) {
		for (Individual individual : individuals) {
			individual.fitness = Fitness.evaluate(individual, this);
		}
	}

	// Sets the fitness of the individuals passed as parameters (whole population)
	private void evaluateIndividual(Individual individual) {
		individual.fitness = Fitness.evaluate(individual, this);
	}


	 //Returns best individual index
	private int getBest() {
		int bestidx = -1;
		for(int i=0; i<population.size() ; i++){
			if (bestidx == -1) {
				bestidx = i;
			} else if (population.get(i).fitness < population.get(bestidx).fitness) {
				bestidx = i;
			}
		}
		return bestidx;
	}

//  -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -

	 // Generates a randomly initialised population
	private ArrayList<Individual> initialise() {
		population = new ArrayList<>();
		for (int i = 0; i < Parameters.popSize; ++i) {
			//chromosome weights are initialised randomly in the constructor
			Individual individual = new Individual();
			population.add(individual);
		}
		evaluateIndividuals(population);
		return population;
	}


	// Generates a randomly initialised population
	private ArrayList<Individual> initialiseAll0() {
		population = new ArrayList<>();
		for (int i = 0; i < Parameters.popSize; ++i) {
			//chromosome weights are initialised randomly in the constructor
			Individual individual = new Individual();
			// predictably set the genes
			for(int gene=0; gene<individual.chromosome.length; gene++){
				individual.chromosome[gene] = 0;
				//individual.chromosome[gene] = (((Parameters.minGene*(-1)) + Parameters.maxGene) / (gene+1)) - (Parameters.minGene*(-1));
			}
			population.add(individual);
		}
		evaluateIndividuals(population);
		return population;
	}

//  -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -

	// Selection
	private ArrayList<Individual> select() {
		Individual parent1 = population.get(Parameters.random.nextInt(Parameters.popSize));
		Individual parent2 = population.get(Parameters.random.nextInt(Parameters.popSize));
		ArrayList parents = new ArrayList();
		parents.add(parent1);
		parents.add(parent2);
		return parents;
	}


	private int selectRW() {
			double totalPopF = 0;
			Individual current;
			double partialFSum = 0;
			// get total fitness of pop
			for (int i = 0; i < Parameters.popSize; i++) {
				current = population.get(i);
				totalPopF += current.fitness;
			}
			double rand = (Parameters.random.nextDouble() * totalPopF);
			// get probability of each individual
			for (int i = 0; i < Parameters.popSize; i++) {
				current = population.get(i);
				partialFSum = partialFSum + current.fitness;
				// select individual
				if (partialFSum >= rand) {
					return i;
				}
			}
		return Parameters.random.nextInt(Parameters.popSize);  // rando
		//return getBest();  // idx of best
	}


	private Individual selectBest() {
		Individual secondBest = population.get(getWorstIndex()).copy();
		Individual best = population.get(getWorstIndex()).copy();
		double shitFit = 100;
		for (int i=0; i<population.size(); i++) {
			if (population.get(i).fitness < shitFit) {
				secondBest = best.copy();
				best = population.get(i).copy();
				shitFit = best.fitness;
			}
		}
		best = population.get(getBest()).copy();
		Individual parent = best.copy();
		if(Parameters.random.nextBoolean()) {
			parent = secondBest.copy();
		}
		return parent.copy();
	}


	private Individual selectWorst() {
		Individual secondWorst = population.get(getBest()).copy();
		Individual worst = population.get(getBest()).copy();
		double worstFit = 0;
		for (int i=0; i<population.size(); i++) {
			if (population.get(i).fitness > worstFit) {
				secondWorst = worst.copy();
				worst = population.get(i).copy();
				worstFit = worst.fitness;
			}
		}
		worst = population.get(getWorstIndex()).copy();
		Individual parent = worst.copy();
		if(Parameters.random.nextBoolean()) {
			parent = secondWorst.copy();
		}
		return parent.copy();
	}

	//tournament w/&w/o RW, and w/&w/o dups' in TPop
	private int selectTournament() {
		int popBestIdx = -1; // index in pop of best
		int TSize = 17;
		ArrayList TPopIdx = new ArrayList(); // list of the indexes of chosen pop inds
		for(int i=0; i<TSize; i++) {   // add randos
			int idx = Parameters.random.nextInt(Parameters.popSize);
			while(TPopIdx.contains(idx)){ // dup check
				idx = Parameters.random.nextInt(Parameters.popSize);
			}
			TPopIdx.add(idx);
		}
		// SELECT BEST  ///
//		for(int i=0; i<TPopIdx.size() ; i++){  // get the best
//			if(popBestIdx==-1){
//				popBestIdx = (int)TPopIdx.get(i); //set to first object in TPopIdx (the index in pop)
//			}
//			else if(population.get((int)TPopIdx.get(i)).fitness < population.get(popBestIdx).fitness){
//				popBestIdx = (int)TPopIdx.get(i);
//			}
//		}
		// ROULETTE WHEEL TOURNAMENT
		double totalPopF = 0;
		int currentIdx = -1;
		double partialFSum = 0;
		// get total fitness of pop
		for (int i = 0; i < TPopIdx.size(); i++) {
			currentIdx = (int)TPopIdx.get(i);
			totalPopF += population.get(currentIdx).fitness;
		}
		double rand = (Parameters.random.nextDouble() * totalPopF);
		// get probability of each individual
		for (int i = 0; i < TPopIdx.size(); i++) {//
			partialFSum += population.get((int)TPopIdx.get(i)).fitness;
			// select individual
			if (partialFSum >= rand) {//
				return (int)TPopIdx.get(i);
			}
		}
		return (int)TPopIdx.get(0);
	}

//  -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -

	/*** Crossover / Reproduction */
	private ArrayList<Individual> crossover(Individual parent1, Individual parent2) {
		ArrayList<Individual> children = new ArrayList<>();
		children.add(parent1.copy());
		children.add(parent2.copy());
		return children;
	}


	private ArrayList<Individual> crossover1Point(Individual Parent1, Individual Parent2){
	    ArrayList<Individual> children = new ArrayList<>();
        int childrenAmount = 2;

        for (int i=0; i<childrenAmount; i++) {
            double split = (Parameters.random.nextDouble() * ((Parent1.chromosome.length) + 1));
            Individual child = new Individual();
            for (int j=0; j<Parent1.chromosome.length; j++){
                if (j<split) {
                    child.chromosome[j] = Parent1.chromosome[j];
                }
                else{
                    child.chromosome[j] = Parent2.chromosome[j];
                }
            }
            children.add(child);
        }
    	return children;
    }


	private ArrayList<Individual> crossover2Point(Individual Parent1, Individual Parent2){
		ArrayList<Individual> children = new ArrayList<>();
		int childrenAmount = 2;
		for (int i=0; i<childrenAmount; i++) {
			double split1 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length)));
			double split2 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length) + 1));
			while ( split2<split1){
				split2 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length) + 1));
			}
			Individual child = new Individual();
			for (int j=0; j<Parent1.chromosome.length; j++){
				if (j<split1) {
					child.chromosome[j] = Parent1.chromosome[j];
				}
				else if (j<split2){
					child.chromosome[j] = Parent2.chromosome[j];
				}
				else{
					child.chromosome[j] = Parent1.chromosome[j];
				}
			}
			children.add(child);
		}
		return children;
	}


    private ArrayList<Individual> crossover3Point(Individual Parent1, Individual Parent2){
        ArrayList<Individual> children = new ArrayList<>();
        int childrenAmount = 2;
        for (int i=0; i<childrenAmount; i++) {
            double split1 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length) - 1));
            double split2 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length)));
            double split3 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length) + 1));
            while ( split2<split1){
                split2 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length)));
            }
            while ( split3<split2){
                split3 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length) + 1));
            }
            Individual child = new Individual();
            for (int j=0; j<Parent1.chromosome.length; j++){
                if (j<split1) {
                    child.chromosome[j] = Parent1.chromosome[j];
                }
                else if (j<split2){
                    child.chromosome[j] = Parent2.chromosome[j];
                }
                else if (j<split3){
                    child.chromosome[j] = Parent1.chromosome[j];
                }
                else{
                    child.chromosome[j] = Parent2.chromosome[j];
                }
            }
            children.add(child);
        }
        return children;
    }


	private ArrayList<Individual> crossover3PointRandom(Individual Parent1, Individual Parent2){
		ArrayList<Individual> children = new ArrayList<>();
		int childrenAmount = 2;
		for (int i=0; i<childrenAmount; i++) {
			double split1 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length) - 1));
			double split2 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length)));
			double split3 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length) + 1));
			while ( split2<split1){
				split2 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length)));
			}
			while ( split3<split2){
				split3 = (Parameters.random.nextDouble() * ((Parent1.chromosome.length) + 1));
			}
			boolean rand1 = Parameters.random.nextBoolean();
			boolean rand2 = Parameters.random.nextBoolean();
			boolean rand3 = Parameters.random.nextBoolean();
			boolean rand4 = Parameters.random.nextBoolean();
			Individual child = new Individual();
			for (int j=0; j<Parent1.chromosome.length; j++){
				if (j<split1) {
					if(rand1) {
						child.chromosome[j] = Parent1.chromosome[j];
					}
					else {
						child.chromosome[j] = Parent2.chromosome[j];
					}
				}
				else if (j<split2){
					if(rand2) {
						child.chromosome[j] = Parent1.chromosome[j];
					}
					else {
						child.chromosome[j] = Parent2.chromosome[j];
					}
				}
				else if (j<split3){
					if(rand3) {
						child.chromosome[j] = Parent1.chromosome[j];
					}
					else {
						child.chromosome[j] = Parent2.chromosome[j];
					}
				}
				else{
					if(rand4) {
						child.chromosome[j] = Parent1.chromosome[j];
					}
					else {
						child.chromosome[j] = Parent2.chromosome[j];
					}
				}
			}
			children.add(child);
		}
		return children;
	}


	private ArrayList<Individual> crossoverLayer(Individual Parent1, Individual Parent2){
		ArrayList<Individual> children = new ArrayList<>();
		int childrenAmount = 2;
		for (int i=0; i<childrenAmount; i++) {
			double split1 = 5 * Parameters.getNumHidden() + Parameters.getNumHidden();
			Individual child = new Individual();
			for (int j=0; j<Parent1.chromosome.length; j++){
				if (j<split1) {
					child.chromosome[j] = Parent1.chromosome[j];
				}
				else {
					child.chromosome[j] = Parent2.chromosome[j];
				}
			}
			children.add(child);
		}
		return children;
	}


    private ArrayList<Individual> crossoverUniform(Individual Parent1, Individual Parent2) {
        ArrayList<Individual> children = new ArrayList<>();
        int childrenAmount = 2; //Math.toIntExact((Parameters.popSize-2)/10);
        for (int i = 0; i < childrenAmount; i++) {
        Individual child = new Individual();
            for (int j=0; j<Parent1.chromosome.length; j++) {
                if(Parameters.random.nextBoolean()){
                    child.chromosome[j] = Parent1.chromosome[j];
                }
                else{
                    child.chromosome[j] = Parent2.chromosome[j];
                }
            }
                children.add(child);
        }
        return children;
    }

//  -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -

	/** Mutation  */
	private void mutate(ArrayList<Individual> individuals) {
		for(Individual individual : individuals) {
			for (int i = 0; i < individual.chromosome.length; i++) {
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {
						individual.chromosome[i] += (Parameters.mutateChange);
					} else {
						individual.chromosome[i] -= (Parameters.mutateChange);
					}
				}
			}
		}
	}

//	private void mutateFlip(ArrayList<Individual> individuals) {
//		for (Individual individual : individuals) {
//			for (int i = 0; i < individual.chromosome.length; i++) {
//				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
//					individual.chromosome[i] = individual.chromosome[i] * (-1);
//				}
//			}
//		}
//	}
//	private void mutateSwap(ArrayList<Individual> individuals) {
//		for (Individual individual : individuals) {
//			int first = (int)(Math.round(Parameters.random.nextDouble()*Parameters.popSize));
//			int second = (int)(Math.round(Parameters.random.nextDouble()*Parameters.popSize));
//			double secondV = individual.chromosome[second];
//			individual.chromosome[second] = individual.chromosome[first];
//			individual.chromosome[first] = secondV;
//		}
//	}

	private void mutateUntilBetter(ArrayList<Individual> individuals) {
		for(Individual individual : individuals) {
			boolean better = false;
			while (!better){
				evaluateIndividual(individual);
				double FBefore = individual.fitness;
				for (int i = 0; i < individual.chromosome.length; i++) {
					if (Parameters.random.nextDouble() < Parameters.mutateRate) {
						if (Parameters.random.nextBoolean()) {
							individual.chromosome[i] += (Parameters.mutateChange);
						} else {
							individual.chromosome[i] -= (Parameters.mutateChange);
						}
					}
				}
				evaluateIndividual(individual);
				if(individual.fitness < FBefore){ //population.get(getBest()).fitness){
					better = true;
				}
			}
		}
	}


	private void mutatePerPair(ArrayList<Individual> individuals) {
		for(Individual individual : individuals) {
			for (int i = 0; i < individual.chromosome.length; i=i+2) {
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {
						individual.chromosome[i] += (Parameters.mutateChange);
						if(i+1<individual.chromosome.length) {
							individual.chromosome[i + 1] += (Parameters.mutateChange);
						}
					} else {
						individual.chromosome[i] -= (Parameters.mutateChange);
						if(i+1<individual.chromosome.length) {
							individual.chromosome[i + 1] -= (Parameters.mutateChange);
						}
					}
				}
			}
		}
	}

	private void mutateByLayer(ArrayList<Individual> individuals) {
		for(Individual individual : individuals) {
			for (int i = 0; i < (Parameters.getNumHidden()); i++) { // first 5H bit
				int idx = 0;
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {
						for(int j=0; j<Parameters.getNumHidden(); j++) {
							individual.chromosome[i+idx] += (Parameters.mutateChange);
							idx=idx+5;
						}
					} else {
						for(int j=0; j<Parameters.getNumHidden(); j++) {
							individual.chromosome[i+idx] -= (Parameters.mutateChange);
							idx=idx+5;
						}
					}
				}
			}
			for (int i = 5*Parameters.getNumHidden(); i < (5*Parameters.getNumHidden()+3); i++) { // hidden bias bit
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {
						individual.chromosome[i] += (Parameters.mutateChange);
					} else {
						individual.chromosome[i] -= (Parameters.mutateChange);
					}
				}
			}
			int idx = 5*Parameters.getNumHidden()+Parameters.getNumHidden();
			for (int i = 0; i < Parameters.getNumHidden(); i++) { // second 3H bit
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {  // +
						individual.chromosome[i+idx] += (Parameters.mutateChange);
						individual.chromosome[i+idx+Parameters.getNumHidden()] += (Parameters.mutateChange);
						individual.chromosome[i+idx+2*Parameters.getNumHidden()] += (Parameters.mutateChange);
						}
					} else {  // -
						individual.chromosome[i+idx] += (Parameters.mutateChange);
						individual.chromosome[i+idx+Parameters.getNumHidden()] += (Parameters.mutateChange);
						individual.chromosome[i+idx+2*Parameters.getNumHidden()] += (Parameters.mutateChange);
					}
				}
			for (int i=(5*Parameters.getNumHidden())+Parameters.getNumHidden()+(3*Parameters.getNumHidden()); i<Parameters.getNumGenes(); i++) { // output bias bit
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {
						individual.chromosome[i] += (Parameters.mutateChange);
					} else {
						individual.chromosome[i] -= (Parameters.mutateChange);
					}
				}
			}
		}
	}

//  -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -

	//replace the worst guys
	private void replace(ArrayList<Individual> individuals) {
		for(Individual individual : individuals) {  // for each kid
			int idx = getWorstIndex();   // get worst ind
			population.set(idx, individual);   // replace idx with ind
		}
	}


	//replace the worst guy only if the kid is better
	private void replaceIfBetter(ArrayList<Individual> individuals) {
		for(Individual individual : individuals) {  // for each kid
			double worst = population.get(getWorstIndex()).fitness;
			if(individual.fitness < worst) {
				population.set(getWorstIndex(), individual);   // replace idx with ind
			}
		}
	}


	// replace all the guys with new randos, except the kid, and the parents
	private void replaceAllWithNewIfBest(ArrayList<Individual> individuals) {
		ArrayList<Individual> newPop = initialiseAll0();  // new pop of 0's
		int i = 0;
		for(Individual individual : individuals) { // add the kids
			//if(individual.fitness > population.get(getBest()).fitness){  // if kids not best then skip it
			//	continue;
			//}
			newPop.set(i, individual);
			i++;
		}
		evaluateIndividuals(newPop);
		population = newPop;
	}


	// replace all with kid clones if the kid is the best
	private void replaceAllWithGoldenChild(ArrayList<Individual> individuals) {
		ArrayList<Individual> newPop = new ArrayList<>();  // new pop
		for(Individual individual : individuals) { // add the kids
			if(individual.fitness > population.get(getBest()).fitness){  // if kids not best then skip it
				continue;
			}
			for (int j = 0; j < Parameters.popSize; j++) {
				newPop.add(individual);  // clone god child
			}
			evaluateIndividuals(newPop);
			population = newPop;  // replace all with god child clones
		}
	}

//  -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -    -  -

	/** Returns the index of the worst member of the population
	 * @return */
	private int getWorstIndex() {
		Individual worst = null;
		int idx = -1;
		for (int i = 0; i < population.size(); i++) {
			Individual individual = population.get(i);
			if (worst == null) {
				worst = individual;
				idx = i;
			} else if (individual.fitness > worst.fitness) {
				worst = individual;
				idx = i; 
			}
		}
		return idx;
	}	

	@Override
	public double activationFunction(double x) {
		if (x < -20.0) {
			return -1.0;
		} else if (x > 20.0) {
			return 1.0;
		}
		return Math.tanh(x);
	}
}
