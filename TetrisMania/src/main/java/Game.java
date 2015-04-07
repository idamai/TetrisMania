import javax.security.auth.callback.Callback;

public class Game {

  public State state;
  protected boolean useGui;
  protected TFrame screen;

  protected long delay;


  public interface Callback {
    public int[] execute(Game g, State s);
  }


  public Game(boolean useGui, long delay) {
    init(useGui, delay);
  }


  private void init(boolean useGui, long delay) {
    this.useGui = useGui;
    this.delay = delay;
    state = new State();
    screen = new TFrame(state);
  }


  private void drawFrame() {
    if (useGui) {
      state.draw();
      state.drawNext(0, 0);
    }
  }


  public int run(Callback callback) {
    while (!gameOver()) {
      int[] nextMove = callback.execute(this,state);
      tick(nextMove);
      Core.sleep(delay);
    }
    return score();
  }


  public void tick(int[] nextMove) {
    state.makeMove(nextMove);
  }


  public boolean gameOver() {
    return state.hasLost();
  }


  public int score() {
    return state.getRowsCleared();
  }

}
