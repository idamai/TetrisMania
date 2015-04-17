import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class MasterSlaveTest {
  Master master;
  Slave slave1;

  @Before
  public void setUp() throws Exception {
    master = new Master(9000);
    slave1 = new Slave("localhost", 9000, 9001);
    master.start();
    slave1.start();
  }

  @After
  public void tearDown() throws Exception {
    master.stop();
    slave1.stop();
  }

  @Test
  public void testProcess() throws Exception {
    Core.sleep(5000);
  }
}