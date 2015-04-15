import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

/**
 *
 * Player skeleton
 */
public class PlayerSkeleton {
	private static PlayerSkeleton _instance = null;
	private static int NUM_OF_RANDOM_CHROMOSOME = 10000;
	private static int LOCAL_MAXIMA_THRESHOLD = 5;
	private static double ACCEPTABLE_SCORE_COEFF = 0.8;
	private static double SAMPLING_COEFFICIENT = 0.1;
	private static Random RANDOM_GENERATOR = new Random();
	private static double MIN_MUTATION_IDX = 0.8;
	private static int highscoreRowCleared;
	private static String LOG_ROWS_CLEARED = "Highscore: %1$s. Turn: %2$s";

	// private double[] weights = new double[22];

	private final static Logger LOGGER = Logger.getLogger(PlayerSkeleton.class
			.getName());

	private PlayerSkeleton() {
	}

	public static PlayerSkeleton getInstance() {
		if (_instance == null) {
			_instance = new PlayerSkeleton();
			return _instance;
		} else {
			return _instance;
		}
	}

	/**
	 * Agent's Strategy: picks a move (horizontal positioning and rotation
	 * applied to the falling object) that maximises (reward + heuristic
	 * function)
	 *
	 * @param s
	 * @param legalMoves
	 *            int[n][a] where n is number of possible moves and a is the
	 *            index of orient and action
	 *
	 * @return move int array [orient, slot]
	 *
	 */
	public int[] pickMove(State s, double[] weights, int[][] legalMoves) {
		// indices for legalMoves from State class
		final int ORIENT = State.ORIENT;
		final int SLOT = State.SLOT;

		// for loop variables
		int orient = -1;
		int slot = -1;
		int[] currentMove = new int[2];
		int[] bestMove = new int[2];
		int nextMoveRowsCleared = 0;
		int currentReward = 0;
		double currentHeuristic = 0;
		double currentUtility = 0;
		double bestUtility = Integer.MIN_VALUE;

		// initialises default move to be the first legal move
		CloneState cState = new CloneState(s);
		int[] defaultMove = legalMoves[0];
		cState.tryMakeMove(defaultMove[ORIENT], defaultMove[SLOT]);
		int defaultRowsCleared = cState.getCCleared();

		boolean lastMove = true;
		boolean bestFound = false;

		// Calculate utility for every legal move in the given array
		for (int n = 0; n < legalMoves.length; n++) {
			// Setting variables for calculating utility
			cState = new CloneState(s);
			currentMove = legalMoves[n];
			orient = currentMove[ORIENT];
			slot = currentMove[SLOT];

			// Given the set of moves, try a move which does not make us lose
			if (cState.tryMakeMove(orient, slot)) {
				lastMove = false;
				currentReward = cState.getCCleared();
				currentHeuristic = calculateHeuristic(weights, cState);
				currentUtility = currentHeuristic + currentReward;

				// Keeping track of best utility and the respective action
				if (currentUtility > bestUtility) {
					bestFound = true;
					bestMove[ORIENT] = orient;
					bestMove[SLOT] = slot;
					bestUtility = currentUtility;
					// Records the rows cleared if the current best move is made
					nextMoveRowsCleared = currentReward;
				}
			}
		}

		// Pick the default move if all legal moves cause us to lose or
		// all moves have the same utility values
		if (lastMove || !bestFound) {
			return defaultMove;
		} else {
			return bestMove;
		}
	}

//	/**
//	 * @param nextRowsCleared
//	 * @param nextTurn
//	 */
//	private void updateHighscore(int nextRowsCleared, int nextTurn) {
//
//	}

	// Genetic algorithm
	// Generate Weight Chromosome
	// This function will generate random weights for each features depending
	// on the number of heur functions will be used. Should only be run during
	// initialization process
	// param: number of functions in the problem set
	public static Vector<double[]> generateWeightChromosome(int N) {
		Vector<double[]> generatedWeights = new Vector<double[]>();
		for (int i = 0; i < NUM_OF_RANDOM_CHROMOSOME; i++) {
			double[] weightArray = new double[N];
			for (int j = 0; j < N; j++)
				// all weights should be negativeS
				weightArray[j] = -1 * RANDOM_GENERATOR.nextDouble();
			generatedWeights.add(weightArray);
		}
		return generatedWeights;
	}

	/*
	 * Reproduce function Generate random cutoff point And marry parent x to
	 * parent y
	 *
	 * @param x Weights of first pair
	 *
	 * @param y Weights of second pair
	 *
	 * @param n length of x and y
	 */
	public static double[] Reproduce(double[] x, double[] y, int n) {
		double[] child = new double[n];
		int cutoff = (int) Math.floor(n * RANDOM_GENERATOR.nextDouble());
		// Make sure it won't go out of bounds
		if (cutoff > n - 1)
			cutoff--;
		for (int i = 0; i < cutoff; i++) {
			child[i] = x[i];
		}
		for (int j = cutoff; j < n; j++) {
			child[j] = y[j];
		}
		// now we mutate on of the elements
		double mutationIndex = RANDOM_GENERATOR.nextDouble();
		if (mutationIndex > MIN_MUTATION_IDX) {
			int mutationPosition = (int) Math.floor(n
					* RANDOM_GENERATOR.nextDouble());
			if (mutationPosition > n - 1)
				mutationPosition--;
			double mutatedValue = -1 * RANDOM_GENERATOR.nextDouble();
			child[mutationPosition] = mutatedValue;
		}
		return child;
	}

	// Genetic Algorithm
	// GeneticAlgorithm Function
	// param: chromosome of weights
	public static double[] GeneticAlgorithm(Vector<double[]> weightChromosomes)
			throws IOException {
		PrintWriter writer = new PrintWriter(new FileWriter(
				"geneticalgolog.txt", true), true);
		// fitness function is the test of the game
		// run
		Vector<WeightsFitnessPair> weightChromosomePopulation = new Vector<WeightsFitnessPair>();

		// round fittest will denote current checked round fitness
		int roundFittest = Integer.MIN_VALUE;
		WeightsFitnessPair currentFittestPair = null;

		ExecutorService taskExecutor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		List<Future<WeightsFitnessPair>> results;
		final Vector<double []> finalWeightChromosomes = weightChromosomes;
		List<Callable<WeightsFitnessPair>> tasks = new ArrayList<Callable<WeightsFitnessPair>>(weightChromosomePopulation.size());
		// initialize population and round fittest
		for (int i = 0; i < finalWeightChromosomes.size(); i++) {
			final int finalI = i;

			tasks.add(new Callable<WeightsFitnessPair>() {
				@Override
				public WeightsFitnessPair call() {
					int weightsFitness = runState(finalWeightChromosomes.get(finalI));
					return new WeightsFitnessPair(
							finalWeightChromosomes.get(finalI), weightsFitness);
				}
			});
		}
		try {
			System.out.println("Tasks ran initial: " + tasks.size());
			results = taskExecutor.invokeAll(tasks, Long.MAX_VALUE, TimeUnit.SECONDS);
			System.out.println("Tasks ran after initial: " + tasks.size());
			System.out.println("Results size after initial: " + results.size());
			for (Future<WeightsFitnessPair> it : results){
				if (currentFittestPair == null ||  it.get().getFitness() > roundFittest){
					currentFittestPair = it.get();
					roundFittest = it.get().getFitness();
				}
				weightChromosomePopulation.add(it.get());
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
		writer.println("Original Fittest:");
		writer.println(currentFittestPair.toString());
		writer.flush();
		int localMaximaRetry = 0;
		// Re run until found the fittest
		// try to escape from local maxima by allowing degrading by going down
		// 10%
		int counter = 0;
		final WeightsFitnessPair finalCurrentFittestPair = currentFittestPair;
		WeightsFitnessPair roundFittestPair = currentFittestPair;
		while (roundFittest >= (ACCEPTABLE_SCORE_COEFF * currentFittestPair.getFitness())
				&& localMaximaRetry < LOCAL_MAXIMA_THRESHOLD) {
			Vector<WeightsFitnessPair> newPopulation = new Vector<WeightsFitnessPair>();
			roundFittest = Integer.MIN_VALUE;
			//parallelize stuffs to run in each round
			final Vector<WeightsFitnessPair> finalWeightChromosomePopulation = weightChromosomePopulation;

			tasks = new ArrayList<Callable<WeightsFitnessPair>>(weightChromosomePopulation.size());
			for (int i = 0; i < finalWeightChromosomePopulation.size(); i++) {
				tasks.add(new Callable<WeightsFitnessPair>() {
					@Override
					public WeightsFitnessPair call() {
						double[] childWeights = ProduceChild(finalCurrentFittestPair,
								finalWeightChromosomePopulation);
						int childFitness = runState(childWeights);
						return new WeightsFitnessPair(childWeights, childFitness);
					}
				});
			}

			try {
				results = taskExecutor.invokeAll(tasks, Long.MAX_VALUE, TimeUnit.SECONDS);
				System.out.println("Tasks ran: " + tasks.size());
				System.out.println("Results size: " + results.size());
				for (Future<WeightsFitnessPair> future: results){
					if (roundFittestPair == null || future.get().getFitness() > roundFittestPair.getFitness()) {
						roundFittestPair = future.get();
					}
					newPopulation.add(future.get());
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
			if (roundFittestPair != null && roundFittestPair.getFitness() > currentFittestPair.getFitness()) {
				currentFittestPair = roundFittestPair;
				localMaximaRetry = 0;
			} else
				localMaximaRetry++;
			counter++;
			weightChromosomePopulation = newPopulation;
			writer.println("Fittest after round " + counter + ":");
			writer.println(currentFittestPair.toString());
			writer.flush();
		}
		writer.println("Exiting GA");
		writer.close();
		return currentFittestPair.getWeights();
	}

	// Tournament runner, returns child
	// Limitation = original population must be bigger than 50
	private static double[] ProduceChild(WeightsFitnessPair currentFittestPair,
			Vector<WeightsFitnessPair> currentPopulation) {
		// Select 10% of population and do a tournament of 2 fittest individual
		int numToBeTaken = (int) Math.ceil(SAMPLING_COEFFICIENT * currentPopulation.size());
		Vector<WeightsFitnessPair> tournamentSample = new Vector<WeightsFitnessPair>();
		for (int i = 0; i < numToBeTaken; i++) {
			int candidateIdx = (int) Math.floor(currentPopulation.size()
					* RANDOM_GENERATOR.nextDouble());
			if (candidateIdx > (currentPopulation.size() - 1))
				candidateIdx--;
			tournamentSample.add(currentPopulation.get(candidateIdx));
		}
		// probability to use best weights directly
		double probIdx = RANDOM_GENERATOR.nextDouble();
		int midPoint = numToBeTaken / 2;
		WeightsFitnessPair parentA;
		if (probIdx > 0.95)
			parentA = currentFittestPair;
		else
			parentA = tournamentSelect(new Vector<WeightsFitnessPair>(
					currentPopulation.subList(0, midPoint + 1)));
		WeightsFitnessPair parentB;
		if (probIdx < 0.05)
			parentB = currentFittestPair;
		else
			parentB = tournamentSelect(new Vector<WeightsFitnessPair>(
					currentPopulation.subList(midPoint + 1,
							currentPopulation.size())));
		double[] child = Reproduce(parentA.getWeights(), parentB.getWeights(),
				parentA.getWeights().length);
		return child;

	}

	private static WeightsFitnessPair tournamentSelect(
			Vector<WeightsFitnessPair> currentPopulation) {
		WeightsFitnessPair pairA;
		WeightsFitnessPair pairB;
		if (currentPopulation.size() == 1)
			return currentPopulation.get(0);
		else if (currentPopulation.size() == 2) {
			pairA = currentPopulation.get(0);
			pairB = currentPopulation.get(1);
		} else {
			int midPoint = currentPopulation.size() / 2;
			pairA = tournamentSelect(new Vector<WeightsFitnessPair>(
					currentPopulation.subList(0, midPoint + 1)));
			pairB = tournamentSelect(new Vector<WeightsFitnessPair>(
					currentPopulation.subList(midPoint + 1,
							currentPopulation.size())));
		}
		if (pairA.getFitness() > pairB.getFitness())
			return pairA;
		else
			return pairB;
	}

	public void calculateFeature(double[] features, CloneState s) {
		int i;
		int maxHeight = 0;
		int[] height = s.getCTop();
		int numRow = s.getCField().length;
		int numCol = s.getCField()[0].length;
		for (i = 0; i < numCol; i++) { // copy the number of row
			features[i] = height[i];
			maxHeight = Math.max(height[i], maxHeight);
			if (i + 1 < numCol)
				features[i + numCol] = Math.abs(height[i + 1] - height[i]);
		}
		features[i++] = maxHeight;
		features[i] = 0;
		for (int j = 0; j < numCol; j++) {
			for (int k = 0; k < height[j]; k++) {
				if (s.getCField()[k][j] == 0) {
					features[i]++;
				}
			}
		}
		return;
	}

	// Based on the default heuristic mentioned in project assignment.
	// features 0 - 9 :10 columns height of the wall.
	// features 10 - 18: absolute difference in height of adjacent walls.
	// feature 19: maximum column height.
	// feature 20: number of holes in the wall.
	public double calculateHeuristic(double[] weights, CloneState s) {
		double sum = 0;
		double features[] = new double[21];
		features[0] = 1;
		calculateFeature(features, s);
		for (int i = 0; i < weights.length; i++) {
			sum += weights[i] * features[i];
		}
		return sum;
	}

	public static int runState(final double[] weights) {
		final PlayerSkeleton p = PlayerSkeleton.getInstance();
		Game g = new Game(false, 1); // create headless game with 20ms tick
										// delay

		g.run(new Game.Callback() {
			/** Implement the below **/
			public int[] execute(Game g, State s) {
				int[] nextMove = p.pickMove(s, weights, s.legalMoves());
				return nextMove;
			}
		});
		LOGGER.info("Score: " + g.score());
		return g.score();
	}

	public static void main(String[] args) {
		// TODO: Add in proper logging library instead of System.out.println
		// int score = runState(null);
		// System.out.println("You have completed "+ score +" rows.");
		// run genetic algo here
		Vector<double[]> weightChromosomes = generateWeightChromosome(21);
		try {
			double[] results = GeneticAlgorithm(weightChromosomes);
			for (int i = 0; i < results.length; i++){
				System.out.print(results[i] + " ");
			}
			System.out.println();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("ended");
		System.exit(0);
		return;
	}

	public static void generateTrainData(String[] args) {
		// generateScores("sim_wgts/sim_wgt_space_100.txt",100);
		return;
	}

	public static int generateScores(String sim_wgts_filename, int numSamples) {
		double[] scores = new double[numSamples];
		try {
			File sim_wgts = new File(sim_wgts_filename);
			if (sim_wgts.exists()) {
				Scanner sc = new Scanner(sim_wgts);
				int sampleNum = 0;
				while (sc.hasNext()) {
					String line = new String(sc.nextLine());
					String[] str_wgt = line.split(" ");

					double[] double_wgt = new double[21];
					int i;
					for (i = 0; i < 21; i++) {
						double_wgt[i] = Double.parseDouble(str_wgt[i]);

					}
					// *** function call to play game here: *** //
					// scores[sampleNum] = runState();
					sampleNum++;
				}
				sc.close();
			}
		} catch (FileNotFoundException fnfe) {
			System.out.println(fnfe.getMessage());
		}

		// export scores
		try (PrintStream output = new PrintStream(new File("scores.txt"));) {
			for (int i = 0; i < scores.length; i++) {
				output.println(scores[i]);
			}
			output.close();
		} catch (FileNotFoundException fnfe) {
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

	public String toString() {
		String str = "This weights pair details are as follows:\n";
		for (int i = 0; i < weights.length; i++) {
			str += weights[i];
			if (i != weights.length - 1)
				str += ", ";
		}
		str += "\n With a score of ";
		str += fitness;
		return str;
	}
}