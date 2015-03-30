import java.util.Random;
import java.util.Vector;


public class PlayerSkeleton {
	private static int NUM_OF_RANDOM_CHROMOSOME = 16;
	private static Random RANDOM_GENERATOR = new Random();
	//implement this function to have a working system
	public int pickMove(State s, int[][] legalMoves) {

		return 0;
	}
	
	// Genetic  algorithm
	// Generate Weight Chromosome
	// This function will generate random weights for each features depending
	// on the number of heur functions will be used. Should only be run during
	// initialization process
	// param:  number of functions in 
	public Vector<double[]> generateWeightChromosome(int N) {
		Vector<double[]>  generatedWeights = new Vector<double[]>();
		for (int i = 0; i < NUM_OF_RANDOM_CHROMOSOME; i++){
			double[] weightArray = new double[N];
			for (int j = 0; j< N; j++)
				weightArray[j] = RANDOM_GENERATOR.nextDouble();
			generatedWeights.add(weightArray);
		}
		return generatedWeights;
	}
	
	// Reproduce function
	// Generate random cutoff point
	// And marry parent x to parent y
	public double[] Reproduce(double[] x, double[] y, int N){
		double[] child = new double[N];
		int cutoff = (int) Math.floor(N * RANDOM_GENERATOR.nextDouble());
		for (int i = 0; i <cutoff; i++){
			child[i] = x[i];
		}
		for (int j = cutoff; j <N; j++){
			child[j] = y[j];
		}
		return child;
	}
	
	
    //TODO: Add in calculation of the features based on the state of the environment
    public double calculateFeature(double[] features){
        return 0.0;
    }

    public double calculateHeuristic(double[] weights, State s){
        double sum = 0;
        double features[] = new double[22];
        features[0] = 1;
        calculateFeature(features);
        for (int i = 0; i < weights.length; i++){
            sum += weights[i]*features[i];
        }
        return sum;
    }

  public static int runState(boolean useGui) {
    State s = new State();
		if (useGui) {
      new TFrame(s);
    }
    PlayerSkeleton p = new PlayerSkeleton();
    while(!s.hasLost()) {
      s.makeMove(p.pickMove(s,s.legalMoves()));
      if (useGui) {
        s.draw();
        s.drawNext(0,0);
      }
      try {
        Thread.sleep(300);
      } catch (InterruptedException e) {
        e.printStackTrace();
      }
    }
    return s.getRowsCleared();
  }

	public static void main(String[] args) {
	    //TODO: Add in proper logging library instead of System.out.println
	    int score = runState(false);
		System.out.println("You have completed "+ score  +" rows.");
		
	}
		
	
}
