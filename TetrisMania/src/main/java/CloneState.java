import java.util.logging.Logger;

/**
 * Created by joel on 4/7/15.
 */
public class CloneState extends State {
	private final static Logger LOGGER = Logger.getLogger(CloneState.class.getName()); 
	// cloned variables for determining utility value
	private int cTurn;
	private int cCleared;
	private int[][] cField;
	private int[] cTop;
	private int cNextPiece;
	private int[][] cLegalMoves;

	CloneState(State original) {
		cTop = original.getTop().clone();
		cField = deepClone2D(original.getField());
		cCleared = original.getRowsCleared();
		cTurn = original.getTurnNumber();
		cNextPiece = original.getNextPiece();
		cLegalMoves = deepClone2D(original.legalMoves());
	}

	private int[][] deepClone2D(int[][] array) {
		int[][] copy = array.clone();
		for (int i = 0; i < copy.length; i++) {
			copy[i] = array[i].clone();
		}
		return copy;
	}

	//returns false if you lose - true otherwise
	public boolean tryMakeMove(int orient, int slot) {
		cTurn++;
		//height if the first column makes contact
		int height = cTop[slot]-getpBottom()[cNextPiece][orient][0];

		//for each column beyond the first in the piece
		for(int c = 1; c < pWidth[cNextPiece][orient];c++) {
			height = Math.max(height,cTop[slot+c]-getpBottom()[cNextPiece][orient][c]);
		}


		//check if game ended
		if(height+getpHeight()[cNextPiece][orient] >= ROWS) {
			lost = true;
			return false;
		}

		//for each column in the piece - fill in the appropriate blocks
		for(int i = 0; i < pWidth[cNextPiece][orient]; i++) {

			//from bottom to cTop of brick
			for(int h = height+getpBottom()[cNextPiece][orient][i]; h < height+getpTop()[cNextPiece][orient][i]; h++) {
				cField[h][i+slot] = cTurn;
			}
		}

		//adjust cTop
		for(int c = 0; c < pWidth[cNextPiece][orient]; c++) {
			cTop[slot+c]=height+getpTop()[cNextPiece][orient][c];
		}

		int rowsCleared = 0;

		//check for full rows - starting at the cTop
		for(int r = height+getpHeight()[cNextPiece][orient]-1; r >= height; r--) {
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

	public int getCPiece() {
		return cNextPiece;
	}

	public int[][] cLegalMoves() {
		return cLegalMoves;
	}
}
