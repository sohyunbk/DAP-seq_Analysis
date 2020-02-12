import java.io.*;  // Import the File class
import java.util.*; // Import the Scanner class to read text files

public class ReadFile {
  public static void main(String[] args) {
    try {
      File myObj = new File("Zea_mays.AGPv4.dna.toplevel.fa");
      Scanner myReader = new Scanner(myObj);
      FileWriter myWriter = new FileWriter("Zea_mays.AGPv4.dna.toplevel_OnlyChr.fa");
      int nSwitch = 0;
      while (myReader.hasNextLine()) {
        String data = myReader.nextLine();
        if (data.startsWith(">")==true){
           if(data.startsWith(">B73V4")) nSwitch = 0;
           else nSwitch = 1;}
        if (nSwitch == 1) myWriter.write(data+"\n");

    }

      myReader.close();
      myWriter.close();
      }
      catch (IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }

    }

  }
