import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

public class Detection {
    
    public static double communityDetection(Graph network, int populationSize, int numberOfGeneration, double r, int elite, double crossoverRate, double mutationRate, int selection){
    	//Generate chromosomes for given population size
    	ArrayList<ArrayList<Integer>> population = null;
    	population = generateChromosomes(network, populationSize);
    	
    	//Use selection, crossover and mutation on the population
    	for(int i = 0; i < numberOfGeneration; i++){
    		
    		//Initiation of a new population
        	ArrayList<ArrayList<Integer>> newPopulation = new ArrayList<ArrayList<Integer>>();
        	
        	//Find groups which aren't connected to each other
        	ArrayList<ArrayList<ArrayList<Integer>>> subsets = findSubsets(population);
        	
        	//Rate each chromosome
        	Map<Integer, Double> fitnessScores = new TreeMap<>();
        	ArrayList<Double> fitnessScoresList = new ArrayList<>();
        	for(int j = 0; j < populationSize; j++){
        		double score = communityScore(subsets.get(j), r, network);
        		fitnessScores.put(j, score);
        		fitnessScoresList.add(score);
        	}
        	
        	//Add best chromosomes to new population
        	ArrayList<Integer> elites = getElites(fitnessScores, elite);
        	for(int j = 0; j < elites.size(); j++) {
        		newPopulation.add(population.get(elites.get(j)));
        	}
        	
        	//Select selection method
        	ArrayList<Integer> selectedChrom = new ArrayList<Integer>();
        	switch(selection) {
        	  case 0:
              	//Select chromosomes using roulette selection
          		selectedChrom = rouletteSelection(fitnessScoresList, elites.size());
        	    break;
        	  case 1:
        	    //Select chromosomes using tournament selection
        		selectedChrom = tournamentSelection(fitnessScoresList, elites.size());
        	    break;
        	  default:
        	    //Select chromosomes using tournament selection
        		selectedChrom = rankSelection(fitnessScores, elites.size());
        	}
        	
    		//If random number is less or equal to crossover rate apply uniform crossover for selected chromosomes and create children
    		for(int j = 0; j < populationSize - elites.size() - 1; j+=2) {
    			Random rand = new Random();
    			if (rand.nextDouble() <= crossoverRate) {
	    			ArrayList<ArrayList<Integer>> children = uniformCrossover(population.get(selectedChrom.get(j)), population.get(selectedChrom.get(j+1)));
	    			
	    			//Add children to new population
	    			newPopulation.add(children.get(0));
	    			newPopulation.add(children.get(1));
    			} else {
    				
    				//add selected chromosomes to new population
    				newPopulation.add(population.get(selectedChrom.get(j)));
	    			newPopulation.add(population.get(selectedChrom.get(j+1)));
    			}
        	}
    		
    		//For some sizes of population after adding elites and children there is one free space in population, so there is a need to add one more chromosome
    		if((populationSize - elites.size()) % 2 == 1) {
    			Random rand = new Random();
    			if (rand.nextDouble() <= crossoverRate) {
	    			ArrayList<ArrayList<Integer>> children = uniformCrossover(population.get(selectedChrom.get(selectedChrom.size()-1)), population.get(selectedChrom.size()-2));
	    			newPopulation.add(children.get(0));
    			} else {
    				newPopulation.add(population.get(selectedChrom.get(selectedChrom.size()-1)));
    			}
    		}
    		//Mutate genes on chromosomes
    		population = mutation(newPopulation, network, mutationRate);
    	}
    	
    	//Find groups which aren't connected to each other in the final population
    	ArrayList<ArrayList<ArrayList<Integer>>> subsets = findSubsets(population);
    	
    	//Rate chromosomes and find the best one
    	ArrayList<Double> fitnessScoresList = new ArrayList<>();
    	int max = 0;
    	for(int j = 0; j < populationSize; j++){
    		double score = communityScore(subsets.get(j), r, network);
    		fitnessScoresList.add(score);
    		if(score > fitnessScoresList.get(max)) {
    			max = j;
    		}
    	}
    	
    	//Print found groups in the network
    	System.out.println("Znalezione grupy w sieci:");
    	for(int j = 0; j < subsets.get(max).size(); j++){
    		System.out.print("[");
    		for(int k = 0; k < subsets.get(max).get(j).size(); k++){
    			//Add one to each vertices so that the vertices have original numbering
        		System.out.print(subsets.get(max).get(j).get(k) + 1 + " ");
        	}
    		System.out.println("]");
    	}
    	System.out.println("Liczba grup: " + subsets.get(max).size());
    	System.out.println("Przystosowanie najlepszego chromosomu: " + fitnessScoresList.get(max) + "\n");
    	return fitnessScoresList.get(max);
    }
    
    
    //Mutate genes in each chromosome in population
    public static ArrayList<ArrayList<Integer>> mutation(ArrayList<ArrayList<Integer>> population, Graph network, double mutationRate){
    	Random rand = new Random();
    	for(int i = 0; i < population.size(); i++) {
    		for(int j = 0; j < population.get(i).size(); j++) {
    			
    			//If random number is lesser or equal mutation rate then change value of gene
        		if(rand.nextDouble() <= mutationRate) {
        			
        			//Add all vertices to list
        			ArrayList<Integer> vertix = new ArrayList<>();
        			for(int k = 0; k < population.get(i).size(); k++) {
        				vertix.add(k);
        			}
        			
        			//Find edge between vertices exists in original network and change value of edge
        			do {
        				int newValue = rand.nextInt(vertix.size());
        				if(network.doesEdgeExist(j, vertix.get(newValue))) {
                			population.get(i).set(j, vertix.get(newValue));
                			break;
        				} else {
        					vertix.remove(newValue);
        				}
        				
        			} while(true);
        		}
        	}
    	}
    	
    	return population;
    }
    
    
    //Apply uniform crossover to two chromosomes
    public static ArrayList<ArrayList<Integer>> uniformCrossover(ArrayList<Integer> parent1, ArrayList<Integer> parent2) {
        Random rand = new Random();
        ArrayList<Integer> child1 = new ArrayList<Integer>();
        ArrayList<Integer> child2 = new ArrayList<Integer>();
        ArrayList<ArrayList<Integer>> children = new ArrayList<ArrayList<Integer>>();
        
        for (int i = 0; i < parent1.size(); i++) {
          if (rand.nextBoolean()) {
        	child1.add(parent1.get(i));
        	child2.add(parent2.get(i));
          } else {
        	child1.add(parent2.get(i));
        	child2.add(parent1.get(i));
          }
        }
        children.add(child1);
        children.add(child2);
        
        return children;
      }
    
    
    //Find best chromosomes in population. Returns position of chromosomes in population
    public static <K, V extends Comparable<? super V>> ArrayList<Integer> getElites(Map<K, V> map, int elite) {
    	
    	//Sort chromosomes by fitness value
        List<Entry<K, V>> list = new ArrayList<>(map.entrySet());
        list.sort(Entry.comparingByValue());
        ArrayList<Integer> elites = new ArrayList<>();
        for (Entry<K, V> entry : list) {
            elites.add((Integer) entry.getKey());
        }
        
        //Take given percent of population as elites
        int numberOfElites = Math.floorDiv(elite * map.size(), 100);
        for(int i = 0; i < map.size() - numberOfElites; i++) {
        	elites.remove(0);
        }
        
        return elites;
    }
    
    
    //Apply roulette selection on population
    public static ArrayList<Integer> rouletteSelection(ArrayList<Double> fitnessScores, int elitesSize) {
    	
        //Calculate the sum of all fitness scores
        double totalFitness = 0;
        for (double fitness : fitnessScores) {
          totalFitness += fitness;
        }
        
        //Initialize the cumulative fitness array
        ArrayList<Double> cumulativeFitness = new ArrayList<>();
        double sum = 0;
        for (double fitness : fitnessScores) {
          sum += fitness / totalFitness;
          cumulativeFitness.add(sum);
        }
        
        //Perform roulette selection
        ArrayList<Integer> selectedIndices = new ArrayList<>();
        Random random = new Random();
        for (int i = 0; i < fitnessScores.size() - elitesSize; i++) {
          double randomValue = random.nextDouble();
          for (int j = 0; j < cumulativeFitness.size(); j++) {
            if (randomValue <= cumulativeFitness.get(j)) {
              selectedIndices.add(j);
              break;
            }
          }
        }
        
        return selectedIndices;
      }
    
    
    //Apply tournament selection on population
    public static ArrayList<Integer> tournamentSelection(ArrayList<Double> population, int elitesSize) {
        
    	int tournamentSize = 2;
        Random random = new Random();
        ArrayList<Integer> tournamentResults = new ArrayList<>();
        
        for(int j = 0; j < population.size() - elitesSize; j++ ) {
            ArrayList<Double> tournament = new ArrayList<>();
            ArrayList<Integer> chromosomeIndexes = new ArrayList<>();
            
	        // Fill tournament arraylist with random individuals from population
	        for (int i = 0; i < tournamentSize; i++) {
	        	int randomChrom = random.nextInt(population.size());
	        	chromosomeIndexes.add(randomChrom);
	            tournament.add(population.get(randomChrom));
	        }
	        
	        // Find the fittest individual in the tournament
	        double fittest = tournament.get(0);
            int bestChrom = chromosomeIndexes.get(0);
	        for (int i = 1; i < tournament.size(); i++) {
	            if (tournament.get(i).compareTo(fittest) > 0) {
	                fittest = tournament.get(i);
	                bestChrom = chromosomeIndexes.get(i);
	            }
	        }
	        tournamentResults.add(bestChrom);        
        }
        return tournamentResults;
    }
    
    
    //Apply rank selection on population
    public static ArrayList<Integer> rankSelection(Map<Integer, Double> map, int elitesSize) {
    	
    	//Sort indexes of chromosomes by their function value
    	map = sortMapByValue(map);
    	
    	//Create chromosome ranking
    	Set<Integer> mapKeys = map.keySet();
    	List<Integer> mapKeysList = new ArrayList<Integer>();
    	mapKeysList.addAll(mapKeys);
    	Map<Integer, Integer> newMap = new TreeMap<>();
    	for(int i = 0; i < map.size(); i++) {
    		newMap.put(i, mapKeysList.get(i));
    	}
    	
    	//Initialize the cumulative fitness array
    	ArrayList<Integer> rankFitness = new ArrayList<>();
    	int cumulativeFitness = 0;
    	for(int i = 0; i < map.size(); i++) {
    		cumulativeFitness += i;
    		rankFitness.add(cumulativeFitness);
    	}
    	
    	//Perform rank selection
    	ArrayList<Integer> selectedIndices = new ArrayList<>();
        Random random = new Random();
        for (int i = 0; i < map.size() - elitesSize; i++) {
          double randomValue = random.nextInt(rankFitness.get(rankFitness.size() - 1));
          for (int j = 0; j < map.size(); j++) {
            if (randomValue <= rankFitness.get(j)) {
              selectedIndices.add(newMap.get(j));
              break;
            }
          }
        }
        return selectedIndices;
    }
    
    
    //Sorting map by value
    public static <K, V extends Comparable<? super V>> Map<K, V> sortMapByValue(Map<K, V> map) {
        List<Entry<K, V>> list = new ArrayList<>(map.entrySet());
        list.sort(Entry.comparingByValue());

        Map<K, V> result = new LinkedHashMap<>();
        for (Entry<K, V> entry : list) {
            result.put(entry.getKey(), entry.getValue());
        }

        return result;
    }
    
    
    //Rate chromosome 
    public static double communityScore(ArrayList<ArrayList<Integer>> subsets, double r, Graph network) {
		double CS = 0;
		
		//For each individual group in chromosome create graph and calculate score
		for(int j = 0; j < subsets.size(); j++) {
			
			//Create graph from group
			Graph subMatrix = new Graph(network.getNumVertices());
			for(int k = 0; k < subsets.get(j).size(); k++) {
				for(int l = 0; l < subsets.get(j).size(); l++) {
					if(network.doesEdgeExist(subsets.get(j).get(k), subsets.get(j).get(l))) {
						subMatrix.addEdge(subsets.get(j).get(k), subsets.get(j).get(l));
					}
				}
			}
			
			//Calculate score
			double M = 0;
			int v = 0;
			double mean = 0;
			for(int k = 0; k < subsets.get(j).size(); k++) {
				mean += Math.pow((double) subMatrix.countEdges(subsets.get(j).get(k)), r) / subsets.get(j).size();
				v += subMatrix.countEdges(subsets.get(j).get(k));
			}
			M += mean / subsets.get(j).size();
			CS += M*v;
		}
		
    	return CS;
    }
    
    
    //For every chromosome in population find groups which aren't connected. Returns arraylists of vertices for every chromosome. Each arraylist is a group which isn't connected to any other group
    public static ArrayList<ArrayList<ArrayList<Integer>>> findSubsets(ArrayList<ArrayList<Integer>> population) {
        ArrayList<ArrayList<ArrayList<Integer>>> result = new ArrayList<>();
        for(int i = 0; i < population.size(); i++){
        	ArrayList<ArrayList<Integer>> subsets = new ArrayList<>();
        	
        	//Make graph representation of chromosome
        	Graph subset = new Graph(population.get(0).size());
        	for(int j = 0; j < population.get(0).size(); j++){
            	subset.addEdge(j, population.get(i).get(j));
            }
        	
        	//Find subgraphs
        	subsets = subset.toSubGraph();
        	result.add(subsets);
        }
        
        return result;
    }
    
    
    //Generate population of chromosomes
    public static ArrayList<ArrayList<Integer>> generateChromosomes(Graph network, int populationSize) {
    	ArrayList<ArrayList<Integer>> population = new ArrayList<ArrayList<Integer>>();
    	int length = network.getNumVertices();
    	for(int i = 0; i < populationSize; i++) {
    		ArrayList<Integer> chromosome = new ArrayList<>();
    		for(int j = 0; j < length; j++) {
        		Random random = new Random();
                int gene = random.nextInt(length);
                
                //Checking if connection between edges exist in original network
                while(!network.doesEdgeExist(j, gene)) {
                	gene = random.nextInt(length);
                }
                chromosome.add(gene);
        	}
    		population.add(chromosome);
    	}
    	
		return population;
    }
    
    
    //Read network from file
    private static Graph readNetworkFromFile(String fileName) throws IOException {
        Graph network;
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            String line;
            line = reader.readLine();
            int numberOfLines = Integer.parseInt(line);
            network = new Graph(numberOfLines);
            while ((line = reader.readLine()) != null) {
                String[] parts = line.split(" ");
                
                //Vertices are numbered from one, so to make things easier subtract one from each number
                network.addEdge(Integer.parseInt(parts[0]) - 1, Integer.parseInt(parts[1]) - 1);
            }
        }
        
        return network;
    }
    

    public static void main(String[] args) throws IOException {
    	
    	//Reading network from file
        Graph network = readNetworkFromFile("src/soc-dolphins.txt");
        
        //Start community detection algorithm
        double result = 0;
        double bestResult = 0;
        double cumulativeResult = 0;
        int numberOfIterations = 10;
        for(int i = 0; i < numberOfIterations ; i += 1) {
            result = communityDetection(network, 200, 30, 1.5, 0, 0.8, 0.02, 0);
            cumulativeResult += result;
            if(result > bestResult) {
            	bestResult = result;
            }
        }
        System.out.println("Œredni wynik: " + cumulativeResult/numberOfIterations);
        System.out.println("Najlepszy wynik: " + bestResult);
        
    }
}