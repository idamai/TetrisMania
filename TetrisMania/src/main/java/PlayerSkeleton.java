
public class PlayerSkeleton {

	//implement this function to have a working system
	public int pickMove(State s, int[][] legalMoves) {

		return 0;
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
