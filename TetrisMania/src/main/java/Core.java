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


  public static String[] tokenize(String str) {
    return str.split("\\s*,\\s*");
  }

  public static double[] concatArray(double[] a, double[] b) {
    double[] result = new double[a.length + b.length];
    for (int i = 0; i < a.length; i++) {
      result[i] = a[i];
    }
    for (int i = 0; i < b.length; i++) {
      result[i + a.length] = b[i];
    }
    return result;
  }


  public static Double[] convertStrArrayToDoubleArr(String[] in) {
    Double[] result = new Double[in.length];
    for (int i = 0; i < in.length; i++) {
      result[i] = Double.parseDouble(in[i]);
    }
    return result;
  }


  public static double[] convertBigToSmall(Double[] in) {
    double[] out = new double[in.length];
    for (int i = 0; i < in.length; i++) {
      out[i] = in[i];
    }
    return out;
  }

}
