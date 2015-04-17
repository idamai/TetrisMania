import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class PacketTest {

  @Before
  public void setUp() throws Exception {

  }

  @After
  public void tearDown() throws Exception {

  }

  @Test
  public void testToAndFromBytesWorks() throws Exception {
    String header = "hello";
    String payload = "world";
    String ip = "localhost";
    int port = 8000;
    Packet p = new Packet(ip, port, header, payload);
    Packet p2 = new Packet(p.toBytes());

    assertTrue(p2.getHeader().equals(header));
    assertTrue(p2.getPayload().equals(payload));
    assertTrue(p2.getPort().equals(port));
    assertTrue(p2.getIp().equals(ip));
  }
}