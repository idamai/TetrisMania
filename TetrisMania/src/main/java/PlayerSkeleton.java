import java.io.*;
import java.net.*;
import java.nio.ByteBuffer;
import java.util.*;
import java.util.concurrent.*;
import java.util.logging.Logger;

/**
 * Player skeleton
 */
public class PlayerSkeleton {
  private static PlayerSkeleton _instance = null;
    private static int NUM_OF_RANDOM_CHROMOSOME = 10000;
//  private static int NUM_OF_RANDOM_CHROMOSOME = 10;
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
   * @param legalMoves int[n][a] where n is number of possible moves and a is the
   *                   index of orient and action
   * @return move int array [orient, slot]
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
    final Vector<double[]> finalWeightChromosomes = weightChromosomes;
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
      for (Future<WeightsFitnessPair> it : results) {
        if (currentFittestPair == null || it.get().getFitness() > currentFittestPair.getFitness()) {
          currentFittestPair = it.get();
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
    while (roundFittestPair.getFitness() >= (ACCEPTABLE_SCORE_COEFF * currentFittestPair.getFitness())
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
        for (Future<WeightsFitnessPair> future : results) {
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
    taskExecutor.shutdownNow(); // terminate all remaining running threads
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
    /** Normal operation **/
    if (args.length == 0) {
      // int score = runState(null);
      // System.out.println("You have completed "+ score +" rows.");
      // run genetic algo here
      Vector<double[]> weightChromosomes = generateWeightChromosome(21);
      try {
        double[] results = GeneticAlgorithm(weightChromosomes);
        for (int i = 0; i < results.length; i++) {
          System.out.print(results[i] + " ");
        }
        System.out.println();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      System.out.println("ended");
      System.exit(0);
    }

    /** Else distributed mode **/
    if (args.length < 2) System.exit(1); // invalid
    String type = args[0];

    if (type.equals("master")) {
      /** Distributed mode / MASTER **/
      int port = Integer.parseInt(args[1]);
      Master master = new Master(port);
      try {
        master.start();
        Core.sleep(5000);
        master.display();
        master.runSlaves();
        while (!master.isDone);
        System.out.println("here");
        System.exit(0);
      } catch (SocketException e) {
        e.printStackTrace();
      }

    } else if (type.equals("slave")) {
      /** Distributed mode / SLAVE **/
      if (args.length < 4) System.exit(1); // invalid
      String masterAddr = args[1];
      int masterPort = Integer.parseInt(args[2]);
      int port = Integer.parseInt(args[3]);
      try {
        while (true) {
          Slave slave = new Slave(masterAddr, masterPort, port);
          slave.start();
          while (!slave.isDone) ; // infinite loop, slave never stops!
          slave.stop();
        }
      } catch (SocketException e) {
        e.printStackTrace();
      }
    }

    System.exit(-1);
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


class Packet {
  private String header;
  private String payload;
  private String ip;
  private int port;
  private String origin; // this is only set and not sent when toBytes() is called


  /**
   * Performs deep clone of object
   *
   * @return
   */
  public Packet clone() {
    Packet p2 = new Packet(this.toBytes());
//    p2.origin = this.origin;
    return p2;
  }

  public String getOrigin() {
    return origin;
  }

  public void setOrigin(String origin) {
    this.origin = origin;
  }


  public int getPort() {
    return port;
  }

  public void setPort(int port) {
    this.port = port;
  }


  public String getIp() {
    return ip;
  }

  public void setIp(String ip) {
    this.ip = ip;
  }


  public String getHeader() {
    return header;
  }

  public void setHeader(String header) {
    this.header = header;
  }

  public String getPayload() {
    return payload;
  }

  public void setPayload(String payload) {
    this.payload = payload;
  }

  public Packet(byte[] in) {
    origin = null;
    int headerLen, payloadLen, ipLen;
    byte[] bHeader, bPayload, bIp;
    ByteBuffer bf = ByteBuffer.wrap(in);

    this.port = bf.getInt();
    headerLen = bf.getInt();
    payloadLen = bf.getInt();
    ipLen = bf.getInt();
    bHeader = new byte[headerLen];
    bPayload = new byte[payloadLen];
    bIp = new byte[ipLen];

    bf.get(bHeader, bf.arrayOffset(), headerLen);
    bf.get(bPayload, bf.arrayOffset(), payloadLen);
    bf.get(bIp, bf.arrayOffset(), ipLen);

    try {
      this.header = new String(bHeader, "UTF-8");
      this.payload = new String(bPayload, "UTF-8");
      this.ip = new String(bIp, "UTF-8");
    } catch (UnsupportedEncodingException e) {
    }
  }


  public Packet(String ip, int port, String header, String payload) {
    this.header = header;
    this.payload = payload;
    this.ip = ip;
    this.port = port;
  }


  public byte[] toBytes() {
    byte[] bHeader, bPayload, bIp;
    ByteBuffer buff;
    int size;
    bHeader = header.getBytes();
    bPayload = payload.getBytes();
    bIp = ip.getBytes();
    // 4 ints
    size = bHeader.length + bPayload.length + bIp.length + 4 + 4 + 4 + 4;
    buff = ByteBuffer.allocate(size);
    buff.putInt(port);
    buff.putInt(bHeader.length).putInt(bPayload.length).putInt(bIp.length);
    buff.put(bHeader).put(bPayload).put(bIp);
    return buff.array();
  }


  @Override
  public String toString() {
    return "<PACKET> " +
            "IP=" + ip + " " +
            "PORT=" + port + " " +
            "Header=" + header + " " +
            "Payload=" + payload;
  }
}


class NodeState {
  public static final int READY = 0;
  public static final int BUSY = 1;
  public static final int DONE = 2;
  public static final int EXECUTE_FINAL_START = 3; // only for master
  public static final int EXECUTE_FINAL_BUSY = 4; // only for master
  public static final int EXECUTE_FINAL_DONE = 5; // only for master
  public static final int STOP = 6;
}


class Slave extends SimpleNode {
  private final long PING_INTERVAL = 2000;
  private volatile long refTime;
  private volatile String masterIp;
  private volatile int masterPort;
  private volatile String id;
  private volatile int state;
  public volatile boolean isDone;

  // the first item is score, remaining are weights
  double[] results = null;

  public Slave(String masterIp, int masterPort, int port) {
    super(port);
    this.isDone = false;
    this.masterIp = masterIp;
    this.masterPort = masterPort;
    refTime = 0;
    state = NodeState.READY;
    id = UUID.randomUUID() + "";
    System.out.println("Created slave=" + id);
  }


  private void doPing() {
    if ((System.currentTimeMillis() - refTime) > PING_INTERVAL) {
      refTime = System.currentTimeMillis();
      String payload = id + "," + port;
      this.send(new Packet(masterIp, masterPort, "PING", payload));
    }
  }


  private String prettyId() {
    return "[" + id + "]";
  }


  Thread computeGamethread = new Thread(new Runnable() {
    @Override
    public void run() {
      System.out.println("COMPUTING");
      Vector<double[]> weightChromosomes = PlayerSkeleton.generateWeightChromosome(21);
      try {
        double[] weights = PlayerSkeleton.GeneticAlgorithm(weightChromosomes);
        double[] score = {PlayerSkeleton.runState(weights)};
        // the first item is score, remaining are weights
        results = Core.concatArray(score, weights);
//        for (int i = 0; i < results.length; i++) {
//          System.out.print(results[i] + " ");
//        }
        System.out.println();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      System.out.println("ENDED COMPUTING");
    }
  });


  private String ArrayToCommaString(double[] in) {
    return Arrays.toString(in).replace("[", "").replace("]", "");
  }


  @Override
  protected void process() {
    // ping just in case IP changes during computation
    doPing();

    // process messages in mailbox
    procMessages();

    // The below is a state machine for normal operation.
    // For behavior that is triggered
    // by messages, see procMessages()
    if (state == NodeState.READY) {

    } else if (state == NodeState.BUSY) {
      if (computeGamethread.getState() == Thread.State.TERMINATED) {
        state = NodeState.DONE;
        // send 3 times to minimize possibility of packet loss
        String payload = id + "," + ArrayToCommaString(results);
        this.send(new Packet(masterIp, masterPort, "DONE", payload));
        this.send(new Packet(masterIp, masterPort, "DONE", payload));
        this.send(new Packet(masterIp, masterPort, "DONE", payload));
      }

    } else if (state == NodeState.DONE) {
      this.isDone = true;
    }
  }

  private void procMessages() {
    while (recvQueue.size() > 0) {
      Packet p = recvQueue.poll();
      String header = p.getHeader();

      if (header.equals("RUN") && (state == NodeState.READY)) {
        System.out.println(prettyId() + " initiated [RUN]");
        computeGamethread.start();
        state = NodeState.BUSY;
      }
    }
  }


  @Override
  public SimpleNode start() throws SocketException {
    super.start();
    state = NodeState.READY;
    return this;
  }


  @Override
  public SimpleNode stop() {
    super.stop();
    state = NodeState.STOP;
    return this;
  }
}


class Master extends SimpleNode {

  protected volatile boolean isDone;
  public volatile int state;
  public Map<String, IpPort> slaves;
  public volatile ConcurrentHashMap<String, Boolean> slaveIds;
  public volatile TreeMap<Double, Double[]> results; // the first is score, second is weights
  private int totalRepliesExpected;

  private class IpPort {
    public String ip;

    public int port;

    public IpPort(String ip, int port) {
      this.ip = ip;
      this.port = port;
    }

    @Override
    public String toString() {
      return "IP=" + ip + " port=" + port;
    }

  }


  public Master(int port) {
    super(port);
    isDone = false;
    state = NodeState.STOP;
    slaves = new ConcurrentHashMap<>();
  }


  @Override
  public SimpleNode start() throws SocketException {
    super.start();
    slaves = new ConcurrentHashMap<>();
    state = NodeState.READY;
    return this;
  }


  @Override
  public SimpleNode stop() {
    super.stop();
    state = NodeState.STOP;
    return this;
  }


  private synchronized void printResults() {
    for (Map.Entry<Double, Double[]> ele : results.entrySet()) {
      Double score = ele.getKey();
      Double[] weights = ele.getValue();
      System.out.println("[ " + score + "] = " + Arrays.toString(weights));
    }
  }


  private Map.Entry<Double, Double[]> bestResult() {
    return results.lastEntry();
  }

  private Thread computeFinalThread;

  @Override
  protected void process() {
    super.process();
    procMessages(); // for message-related behavior
    if (state == NodeState.READY) {

    } else if (state == NodeState.BUSY) {
      if (slavesDone()) {
        state = NodeState.DONE;
      }

    } else if (state == NodeState.DONE) {
      System.out.println("ALL DONE");
      // do computation stuff here
      state = NodeState.EXECUTE_FINAL_START;

    } else if (state == NodeState.EXECUTE_FINAL_START) {
      computeFinalThread = new Thread(new Runnable() {
        @Override
        public void run() {
          System.out.println("--- ALL SLAVES COMPLETED ---");
          System.out.println("RESULTS:");
          printResults();

          System.out.println("BEST RESULT:");
          Map.Entry<Double, Double[]> best = bestResult();
          if (best == null) return;
          Double bestScore = best.getKey();
          Double[] bestWeights = best.getValue();

          System.out.println("Executing final iteration..");
          int finalScore = PlayerSkeleton.runState(Core.convertBigToSmall(bestWeights));

          System.out.println("FINAL RESULT");
          System.out.println("Number of slaves: " + slaveIds.size());
          System.out.println("Final Score: " + finalScore);
          System.out.println("Final weights: " + Arrays.toString(bestWeights));

        }
      });


      computeFinalThread.start();
      state = NodeState.EXECUTE_FINAL_BUSY;

    } else if (state == NodeState.EXECUTE_FINAL_BUSY) {
      if (computeFinalThread.getState() == Thread.State.TERMINATED) {
        state = NodeState.EXECUTE_FINAL_DONE;
        System.out.println("Done!!!");
      }

    } else if (state == NodeState.EXECUTE_FINAL_DONE) {
      // TODO add in print out result here
      isDone = true;
    }
  }


  private boolean slavesDone() {
    if (slaveIds == null) {
      return false;
    }

//    if (slaveIds.size() != totalRepliesExpected) {
//      return false;
//    }

    for (boolean isDone : slaveIds.values()) {
      if (!isDone) return false;
    }

    return true;
  }


  private void procMessages() {
    while (hasNext()) {
      Packet p = next();
      String header = p.getHeader();
      if (header.equals("PING")) {
        // sent as <id>,<port>
        String[] tokens = Core.tokenize(p.getPayload());
        String id = tokens[0];
        String port = tokens[1];
        slaves.put(id, new IpPort(p.getOrigin(), Integer.parseInt(port)));

      } else if (header.equals("DONE")) {
        if (state == NodeState.BUSY) {
          String[] args = Core.tokenize(p.getPayload());
          String slaveId = args[0] + "";
          double score = Double.parseDouble(args[1]);
          Double[] weights = Core.convertStrArrayToDoubleArr(Arrays.copyOfRange(args, 2, args.length - 1));
          slaveIds.put(slaveId, true); // this ID is done
          results.put(score, weights);
        }
      }
    }
  }


  /**
   * Command to start slaves running
   */
  public void runSlaves() {
    if (state != NodeState.READY) {
      return;
    }
    state = NodeState.BUSY;
    results = new TreeMap<>(); // get new results on every itr
    getActiveSlaves(); // update active slaves list
    totalRepliesExpected = slaves.size();
    broadcast((Set<String>) slaveIds.keySet(), new Packet("", 0, "RUN", ""));
  }


  /**
   * Broadcasts to a list of slaves, a packet
   *
   * @param slaveIds
   * @param p
   */
  private void broadcast(Set<String> slaveIds, Packet p) {
    for (String id : slaveIds) {
      Packet p2 = p.clone(); // ensure that the refs changed
      IpPort curr = slaves.get(id);
      if (curr != null) { // for safety
        p2.setIp(curr.ip);
        p2.setPort(curr.port);
        send(p2);
      }
    }
  }


  private synchronized void getActiveSlaves() {
    slaveIds = new ConcurrentHashMap<>();
    Iterator it = slaves.keySet().iterator();
    while (it.hasNext()) {
      slaveIds.put((String) it.next(), false); // append slave ID available
    }
  }


  public synchronized String dump() {
    StringBuilder sb = new StringBuilder();
    Iterator it = slaves.entrySet().iterator();
    sb.append("=== Slaves ===\n");
    while (it.hasNext()) {
      Map.Entry curr = (Map.Entry) it.next();
      sb.append("[" + curr.getKey() + "]  " + curr.getValue() + "\n");
    }
    sb.append("=== ===");
    return sb.toString();
  }


  @Override
  public String toString() {
    return dump();
  }


  public void display() {
    System.out.println(dump());
  }
}


class SimpleNode {
  protected int port;
  protected volatile DatagramSocket incoming, outgoing;
  protected boolean isRunning;
  protected volatile ConcurrentLinkedQueue<Packet> sendQueue, recvQueue;

  public SimpleNode(int port) {
    this.port = port;
    isRunning = false;
  }


  public SimpleNode stop() {
    System.out.println("Stopping node..");
    isRunning = false;
    return this;
  }


  public SimpleNode start() throws SocketException {
    if (isRunning) {
      return this;
    }

    isRunning = true;

    incoming = new DatagramSocket(null);
    incoming.setReuseAddress(true);
    incoming.bind(new InetSocketAddress((port)));

    outgoing = new DatagramSocket();

    sendQueue = new ConcurrentLinkedQueue<>();
    recvQueue = new ConcurrentLinkedQueue<>();

    Runnable inBound = new Runnable() {
      @Override
      public void run() {
        while (isRunning) {
          byte[] buffer = new byte[2048];
          Packet p;
          DatagramPacket datagramPacket = new DatagramPacket(buffer, buffer.length);
          try {
            incoming.receive(datagramPacket);
            p = new Packet(buffer);
            p.setOrigin(datagramPacket.getAddress().getHostAddress());
            recvQueue.offer(p);
          } catch (IOException e) {
            e.printStackTrace();
          }
        }
        incoming.close();
      }
    };


    Runnable outBound = new Runnable() {
      @Override
      public void run() {
        while (isRunning) {
          while (sendQueue.size() > 0) {
            Packet p = sendQueue.poll();
            sendPacket(p);
          }
        }
        outgoing.close();
      }
    };


    Runnable processThread = new Runnable() {
      @Override
      public void run() {
        while (isRunning) {
          process();
        }
      }
    };

    new Thread(inBound).start();
    new Thread(outBound).start();
    new Thread(processThread).start();

    Core.sleep(100);

    System.out.println("Node started.");
    return this;
  }


  /**
   * Override this
   */
  protected void process() {

  }


  public void send(Packet p) {
    sendQueue.add(p);
  }


  public Packet next() {
    return recvQueue.poll();
  }


  public boolean hasNext() {
    return !recvQueue.isEmpty();
  }


  private void sendPacket(Packet p) {
    if (outgoing == null) return;
    byte[] data = p.toBytes();
    try {
      outgoing.send(new DatagramPacket(data,
              data.length,
              makeAddr(p.getIp()),
              p.getPort()));
    } catch (IOException e) {
      e.printStackTrace();
    }
  }


  private InetAddress makeAddr(String addr) throws MalformedURLException, UnknownHostException {
    InetAddress result = null;
    try {
      result = InetAddress.getByName(addr);
    } catch (UnknownHostException e) {
      result = InetAddress.getByName(new URL(addr).getHost());
    }
    return result;
  }
}


/**
 * Created by joel on 4/7/15.
 */
class Core {

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


class Game {

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
    // do not instantiate if headless
    if (useGui) {
      screen = new TFrame(state);
    }
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
