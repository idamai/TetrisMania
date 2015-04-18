import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class SimpleNodeTest {

  private SimpleNode node1, node2, node3;

  @Before
  public void setUp() throws Exception {
    node1 = new SimpleNode(5432);
    node2 = new SimpleNode(5433);
    node3 = new SimpleNode(5434);

    node1.start();
    node2.start();
    node3.start();
  }

  @After
  public void tearDown() throws Exception {
    node1.stop();
    node2.stop();
    node3.stop();
  }

  @Test
  public void testSend() throws Exception {
    node1.send(new Packet("localhost", 5433, "hello", "world"));
    Thread.sleep(1000);
//    System.out.println(node1.recvQueue.size());
//    System.out.println(node2.recvQueue.size());
//    System.out.println(node3.recvQueue.size());
    if (node2.hasNext()) {
      System.out.println(node2.next().toString());
    }
  }
}