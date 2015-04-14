import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;
import java.util.Vector;
import java.util.logging.Logger;


/**
 *
 * Player skeleton
 */
public class PlayerSkeleton {
	private static int NUM_OF_RANDOM_CHROMOSOME = 5000;
	private static Random RANDOM_GENERATOR = new Random();
	// private double[] weights = new double[22];

	private final static Logger LOGGER = Logger.getLogger(PlayerSkeleton.class.getName());
	/**
	 * Agent's Strategy: picks a move (horizontal positioning and rotation applied to the falling object)
	 * that maximises (reward + heuristic function)
	 *
	 * @param s
	 * @param legalMoves	int[n][a] where n is number of possible moves and a is the index of orient and action
	 *
	 * @return move			int array [orient, slot]
	 *
	 */
	public int[] pickMove(State s, double[] weights, int[][] legalMoves) {
		final String LOG_BEST_VALUES = "bU: %1$s cU %2$s cH: %3$s cR: %3$s";
		final String LOG_BEST_MOVE = "BEST move picked: %1$s";
		final String LOG_EQUAL_MOVE = "FIRST move picked: %1$s (Equal utility values)";
		final String LOG_LOSING_MOVE = "LOSING move picked: %1$s";

		// indices for legalMoves from State class
		final int ORIENT = State.ORIENT;
		final int SLOT = State.SLOT;

		// initialise variables for finding arg max (utility)
		CloneState cState;
		int orient = -1;
		int slot = -1;
		int currentReward = 0;
		double currentHeuristic = 0;
		double currentUtility = 0;
		double bestUtility = currentUtility;
		int[] currentAction = new int[2];
		int[] bestMove = new int[2];
		boolean bestFound = false;
		boolean lastMove = true;

		// Calculate utility for every legal move in the given array
		for (int n = 0; n < legalMoves.length; n++) {
			// reset values
			currentReward = 0;
			currentHeuristic = 0;
			currentUtility = 0;

			// Setting variables for calculating utility
			cState = new CloneState(s);
			currentAction = legalMoves[n];
			orient = currentAction[ORIENT];
			slot = currentAction[SLOT];

			assert(cState.getCPiece() == s.getNextPiece());
			assert(cState.cLegalMoves()[n][ORIENT] == orient);
			assert(cState.cLegalMoves()[n][SLOT] == slot);
			// Given the set of moves, try a move which does not make us lose
			if (cState.tryMakeMove(orient, slot)) {
				lastMove = false;
				currentReward = cState.getCCleared();
				currentHeuristic = calculateHeuristic(weights, cState);
				currentUtility = currentReward + currentHeuristic;

				// Keeping track of max utility and the respective action
				if (currentUtility > bestUtility) {
					if (!bestFound) {
						bestFound = true;
					}
					LOGGER.info(String.format(LOG_BEST_VALUES, bestUtility, currentUtility, currentHeuristic, currentReward));
					bestMove[ORIENT] = orient;
					bestMove[SLOT] = slot;
					bestUtility = currentUtility;
				}
			}
		}

		// all of the legal moves end the game
		if (lastMove) {
			int[] losingMove = legalMoves[0];
			LOGGER.fine(String.format(LOG_LOSING_MOVE, Arrays.toString(losingMove)));
			return losingMove;
		}

		if (bestFound) {
			LOGGER.fine(String.format(LOG_BEST_MOVE, Arrays.toString(bestMove)));
			return bestMove;
		} else {
			// all moves have some utility value, pick the first legal move
			int[] firstMove = legalMoves[0];
			LOGGER.fine(String.format(LOG_EQUAL_MOVE, Arrays.toString(firstMove)));
			return firstMove;
		}
	}

	// Genetic  algorithm
	// Generate Weight Chromosome
	// This function will generate random weights for each features depending
	// on the number of heur functions will be used. Should only be run during
	// initialization process
	// param:  number of functions in the problem set
	public static Vector<double[]> generateWeightChromosome(int N) {
		Vector<double[]>  generatedWeights = new Vector<double[]>();
		for (int i = 0; i < NUM_OF_RANDOM_CHROMOSOME; i++){
			double[] weightArray = new double[N];
			for (int j = 0; j< N; j++)
				// all weights should be negativeS
				weightArray[j] = -1 * RANDOM_GENERATOR.nextDouble();
			generatedWeights.add(weightArray);
		}
		return generatedWeights;
	}

	/* Reproduce function
	* Generate random cutoff point
	* And marry parent x to parent y
	* @param x Weights of first pair
	* @param y Weights of second pair
	* @param n length of x and y
	*/
	public static double[] Reproduce(double[] x, double[] y, int n){
		double[] child = new double[n];
		int cutoff = (int) Math.floor(n * RANDOM_GENERATOR.nextDouble());
		//Make sure it won't go out of bounds
		if (cutoff > n-1)
			cutoff--;
		for (int i = 0; i <cutoff; i++){
			child[i] = x[i];
		}
		for (int j = cutoff; j <n; j++){
			child[j] = y[j];
		}
		// now we mutate on of the elements
		double mutationIndex = RANDOM_GENERATOR.nextDouble();
		if (mutationIndex > 0.8){
			int mutationPosition =(int) Math.floor(n * RANDOM_GENERATOR.nextDouble());
			if (mutationPosition > n-1)
				mutationPosition--;
			double mutatedValue = -1 * RANDOM_GENERATOR.nextDouble();
			child[mutationPosition]= mutatedValue;
		}
		return child;
	}

	// Genetic Algorithm
	// GeneticAlgorithm Function
	// param: chromosome of weights
	public static double[] GeneticAlgorithm(Vector<double[]> weightChromosomes) throws IOException{
		PrintWriter writer = new PrintWriter(new FileWriter("geneticalgolog.txt", true),true);

		//fitness function is the test of the game
		//run
		Vector<WeightsFitnessPair>  weightChromosomePopulation = new Vector<WeightsFitnessPair>();
		// round fittest will denote current checked round fitness
		int roundFittest=Integer.MIN_VALUE;

		WeightsFitnessPair currentFittestPair = null;

		//initialize population and round fittest
		for (int i = 0; i< weightChromosomes.size(); i++){
			int weightsFitness = runState(weightChromosomes.get(i));
			WeightsFitnessPair weightsFitnessPair = new WeightsFitnessPair(weightChromosomes.get(i), weightsFitness);
			weightChromosomePopulation.add(weightsFitnessPair);
			if (weightsFitness > roundFittest){
				roundFittest = weightsFitness;
				currentFittestPair = weightsFitnessPair;
			}
			//make sure at least its initialized
			if (currentFittestPair == null)
				currentFittestPair= weightsFitnessPair;
		}
		writer.println("Original Fittest:");
		writer.println(currentFittestPair.toString());
		int localMaximaRetry = 0;
		// Re run until found the fittest
		// try to escape from local maxima by allowing degrading by going down 10%
		int counter = 0;
		while (roundFittest >= ( 0.9 * currentFittestPair.getFitness()) && localMaximaRetry < 100 ){
			Vector<WeightsFitnessPair> newPopulation = new Vector<WeightsFitnessPair>();
			roundFittest = Integer.MIN_VALUE;
			WeightsFitnessPair roundFittestPair = null;
			for (int i = 0; i < weightChromosomePopulation.size(); i++){
				//Choose 4 parents candidate
				int candidateAIdx = (int) Math.floor(weightChromosomePopulation.size() * RANDOM_GENERATOR.nextDouble());
				int candidateBIdx = (int) Math.floor(weightChromosomePopulation.size() * RANDOM_GENERATOR.nextDouble());
				int candidateCIdx = (int) Math.floor(weightChromosomePopulation.size() * RANDOM_GENERATOR.nextDouble());
				int candidateDIdx = (int) Math.floor(weightChromosomePopulation.size() * RANDOM_GENERATOR.nextDouble());
				//Ensure none of them out of bounds
				if (candidateAIdx > weightChromosomePopulation.size()-1)
					candidateAIdx--;
				if (candidateBIdx > weightChromosomePopulation.size()-1)
					candidateBIdx--;
				if (candidateCIdx > weightChromosomePopulation.size()-1)
					candidateCIdx--;
				if (candidateDIdx > weightChromosomePopulation.size()-1)
					candidateDIdx--;
				//simulate probability of fitness function using "tournament" styleS
				WeightsFitnessPair candidateA = weightChromosomePopulation.get(candidateAIdx);
				WeightsFitnessPair candidateB = weightChromosomePopulation.get(candidateBIdx);
				WeightsFitnessPair candidateC = weightChromosomePopulation.get(candidateCIdx);
				WeightsFitnessPair candidateD = weightChromosomePopulation.get(candidateDIdx);

				WeightsFitnessPair parentX;
				WeightsFitnessPair parentY;
				if (candidateA.getFitness() > candidateB.getFitness())
					parentX = candidateA;
				else
					parentX = candidateB;

				if (candidateC.getFitness() > candidateD.getFitness())
					parentY = candidateC;
				else
					parentY = candidateD;

				double[] childWeights = Reproduce(parentX.getWeights(), parentY.getWeights(), parentX.getWeights().length);
				int childFitness= runState(childWeights);
				WeightsFitnessPair childFitnessPair = new WeightsFitnessPair(childWeights, childFitness);
				if (childFitness > roundFittest){
					roundFittest = childFitness;
					roundFittestPair = childFitnessPair;
				}
				if (roundFittestPair == null)
					roundFittestPair = childFitnessPair;
				newPopulation.add(childFitnessPair);
			}
			if (roundFittestPair.getFitness() > currentFittestPair.getFitness()){
				currentFittestPair = roundFittestPair;
				localMaximaRetry = 0;
			} else
				localMaximaRetry++;
			counter++;
			weightChromosomePopulation = newPopulation;
			writer.println("Fittest after round "+counter+":");
			writer.println(currentFittestPair.toString());
		}
		writer.close();
		return currentFittestPair.getWeights();

	}


    public void calculateFeature(double[] features, CloneState s){
        int i;
        int maxHeight = 0;
        int[] height = s.getCTop();
		int numRow = s.getCField().length;
		int numCol = s.getCField()[0].length;
        for (i=0; i < numCol; i++){       //copy the number of row
            features[i] = height[i];
            maxHeight = Math.max(height[i], maxHeight);
            if (i+1 < numCol)
                features[i+numCol] = Math.abs(height[i + 1] - height[i]);
        }
        features[i++] = maxHeight;
        features[i] = 0;
        for (int j = 0; j < numCol; j++){
            for (int k = 0; k < height[j];k++){
                if (s.getCField()[k][j] == 0) {
                    features[i]++;
                }
            }
        }
        return;
    }

	//Based on the default heuristic mentioned in project assignment.
	//features 0 - 9 :10 columns height of the wall.
	//features 10 - 18: absolute difference in height of adjacent walls.
	//feature 20: maximum column height.
	//feature 21: number of holes in the wall.
    public double calculateHeuristic(double[] weights, CloneState s){
        double sum = 0;
        double features[] = new double[22];
        features[0] = 1;
        calculateFeature(features, s);
        for (int i = 0; i < weights.length; i++){
            sum += weights[i]*features[i];
        }
        return sum;
    }


  public static int runState(final double[] weights) {
    final PlayerSkeleton p = new PlayerSkeleton();
    Game g = new Game(false, 1); // create headless game with 20ms tick delay

    g.run(new Game.Callback() {
      /** Implement the below **/
      public int[] execute(Game g, State s) {
        int[] nextMove = p.pickMove(s,weights, s.legalMoves());
        return nextMove;
      }
    });

    return g.score();
  }


	public static void main(String[] args) {
	    //TODO: Add in proper logging library instead of System.out.println
	    //int score = runState(null);
		//System.out.println("You have completed "+ score  +" rows.");
		// run genetic algo here
		Vector<double[]> weightChromosomes = generateWeightChromosome(21);
		try {
			GeneticAlgorithm(weightChromosomes);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}


	public static void generateTrainData(String[] args){
	    //generateScores("sim_wgts/sim_wgt_space_100.txt",100);
	    return;
	  }

	  public static int generateScores(String sim_wgts_filename, int numSamples){
	    double[] scores = new double[numSamples];
	    try{
	      File sim_wgts = new File(sim_wgts_filename);
	      if(sim_wgts.exists()){
	        Scanner sc = new Scanner(sim_wgts);
	        int sampleNum = 0;
	        while(sc.hasNext()){
	          String line = new String(sc.nextLine());
	          String[] str_wgt = line.split(" ");

	          double[] double_wgt = new double[21];
	          int i;
	          for(i=0; i < 21; i++){
	            double_wgt[i] = Double.parseDouble(str_wgt[i]);

	          }
	          // *** function call to play game here: *** //
	          //scores[sampleNum] = runState();
	          sampleNum++;
	        }
	        sc.close();
	      }
	    }catch(FileNotFoundException fnfe){
	      System.out.println(fnfe.getMessage());
	    }

	    // export scores
	    try(PrintStream output = new PrintStream(new File("scores.txt"));){
	      for(int i=0; i<scores.length;i++){
	        output.println(scores[i]);
	      }
	      output.close();
	    }catch(FileNotFoundException fnfe){
	      System.out.println(fnfe.getMessage());
	    }

	    return 0;
	  }

}

class WeightsFitnessPair {
	private double[] weights;
	private int fitness;

	public WeightsFitnessPair(double[] weights, int fitness) {
		this.weights = weights.clone();
		this.fitness = fitness;
	}

	public double[] getWeights() {
		return weights;
	}

	public int getFitness() {
		return fitness;
	}

	public String toString(){
		String str = "This weights pair details are as follows:\n";
		for (int i = 0; i<weights.length;i++){
			str+=weights[i];
			if (i!= weights.length-1)
				str+=", ";
		}
		str+="\n With a score of ";
		str+=fitness;
		return str;
	}
}
