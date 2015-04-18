import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class MasterSlaveTest {
  Master master;
  Slave slave1,slave2;

  @Before
  public void setUp() throws Exception {
    master = new Master(9000);
    slave1 = new Slave("localhost", 9000, 9001);
    slave2 = new Slave("localhost", 9000, 9002);
    slave1.start();
    slave2.start();
    master.start();
  }

  @After
  public void tearDown() throws Exception {
    master.stop();
    slave1.stop();
    slave2.stop();
  }

  @Test
  public void testProcess() throws Exception {
    Core.sleep(4000);
    master.display();
    master.runSlaves();
    Core.sleep(7000);
  }
}