/**
 * Created by joel on 4/7/15.
 */
public class Core {

  public static void sleep(long ms) {
    try {
      Thread.sleep(ms);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }
  }

}