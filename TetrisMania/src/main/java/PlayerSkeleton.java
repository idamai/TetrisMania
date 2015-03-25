
public class PlayerSkeleton {

	//implement this function to have a working system
	public int pickMove(State s, int[][] legalMoves) {

		return 0;
	}


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
	
	public static void main(String[] args) {
    //TODO: Add in proper logging library instead of System.out.println
		State s = new State();
		new TFrame(s);
		PlayerSkeleton p = new PlayerSkeleton();
		while(!s.hasLost()) {
			s.makeMove(p.pickMove(s,s.legalMoves()));
			s.draw();
			s.drawNext(0,0);
			try {
				Thread.sleep(300);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		System.out.println("You have completed "+s.getRowsCleared()+" rows.");
	}
	
}
