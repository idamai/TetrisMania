import java.util.Random;
import java.util.Vector;


public class PlayerSkeleton {
	private static int NUM_OF_RANDOM_CHROMOSOME = 16;
	private static Random RANDOM_GENERATOR = new Random();
	private double[] weights;
	
	/**
	 * Agent's Strategy: picks a move (horizontal positioning and rotation applied to the falling object)
	 * that maximises (reward + heuristic function)
	 * 
	 * @param State
	 * @param legalMoves	int[n][a] where n is number of possible moves and a is the index of orient and action
	 * 
	 * @return move			int array [orient, slot]
	 * 
	 */
	public int[] pickMove(State s, int[][] legalMoves) {
		
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
		int[] bestAction = new int[2];
		boolean bestFound = false;
		
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
			
			// Given the set of moves, try a move which does not make us lose
			if (cState.tryMakeMove(orient, slot)) {
				currentReward = cState.getCCleared(); 
				currentHeuristic = calculateHeuristic(cState);
				currentUtility = currentReward + currentHeuristic;
				
				// Keeping track of max utility and the respective action
				if (currentUtility > bestUtility) {
					if (!bestFound) {
						bestFound = true;
					}
					bestAction[ORIENT] = orient;
					bestAction[SLOT] = slot;
				}
			}
		}
		
		// If best move is not available, it means all legal moves are losing moves.
		// The first legal (losing) move is then returned
		if (bestFound) {
			return bestAction;
		} else {
			int[] losingMove = legalMoves[0];
			return losingMove;
		}
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
	
	
    public void calculateFeature(double[] features, CloneState s){
        int i;
        int maxHeight = 0;
        int[] height = s.getCTop();
        int numCol = s.getCField()[0].length;
        int numRow = s.getCField().length;
        for (i=0; i < numCol; i++){       //copy the number of row
            features[i] = height[i];
            maxHeight = Math.max(height[i], maxHeight);
            if (i+1 < numCol)
                features[i+numCol] = Math.abs(height[i + 1] - height[i]);
        }
        features[i++] = maxHeight;
        features[i] = 0;
        for (int j = 0; j < numRow; j++){
            for (int k = 0; k < height[j];k++){
                if (s.getCField()[j][k] == 0) {
                    features[i]++;
                }
            }
        }
        return;
    }

    public double calculateHeuristic(CloneState s){
        double sum = 0;
        double features[] = new double[22];
        features[0] = 1;
        calculateFeature(features, s);
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
		
	public class CloneState extends State {
		// cloned variables for determining utility value
		private int cTurn;	
		private int cCleared;
		private int[][] cField;
		private int[] cTop; 			
		
		CloneState(State original) {
			cTop = original.getTop().clone();
			cField = original.getField().clone();
			cCleared = original.getRowsCleared();
			cTurn = original.getTurnNumber();
		}
		
		//returns false if you lose - true otherwise
		public boolean tryMakeMove(int orient, int slot) {
			cTurn++;
			//height if the first column makes contact
			int height = cTop[slot]-getpBottom()[nextPiece][orient][0];
			//for each column beyond the first in the piece
			for(int c = 1; c < pWidth[nextPiece][orient];c++) {
				height = Math.max(height,cTop[slot+c]-getpBottom()[nextPiece][orient][c]);
			}
			
			//check if game ended
			if(height+getpHeight()[nextPiece][orient] >= ROWS) {
				lost = true;
				return false;
			}

			//for each column in the piece - fill in the appropriate blocks
			for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
				
				//from bottom to cTop of brick
				for(int h = height+getpBottom()[nextPiece][orient][i]; h < height+getpTop()[nextPiece][orient][i]; h++) {
					cField[h][i+slot] = cTurn;
				}
			}
			
			//adjust cTop
			for(int c = 0; c < pWidth[nextPiece][orient]; c++) {
				cTop[slot+c]=height+getpTop()[nextPiece][orient][c];
			}
			
			int rowsCleared = 0;
			
			//check for full rows - starting at the cTop
			for(int r = height+getpHeight()[nextPiece][orient]-1; r >= height; r--) {
				//check all columns in the row
				boolean full = true;
				for(int c = 0; c < COLS; c++) {
					if(cField[r][c] == 0) {
						full = false;
						break;
					}
				}
				//if the row was full - remove it and slide above stuff down
				if(full) {
					rowsCleared++;
					cCleared++;
					//for each column
					for(int c = 0; c < COLS; c++) {

						//slide down all bricks
						for(int i = r; i < cTop[c]; i++) {
							cField[i][c] = cField[i+1][c];
						}
						//lower the cTop
						cTop[c]--;
						while(cTop[c]>=1 && cField[cTop[c]-1][c]==0)	cTop[c]--;
					}
				}
			}
			return true;
		}
		
		public int getCTurn() {
			return cTurn;
		}
		
		public int getCCleared() {
			return cCleared;
		}
		
		public int[][] getCField() {
			return cField;
		}
		
		public int[] getCTop() {
			return cTop;
		}
	}
	
	
	
}
