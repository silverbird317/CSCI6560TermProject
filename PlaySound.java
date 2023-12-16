import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.SourceDataLine;

public class PlaySound {

	private static int frequency = 44100;
	private static int durationMs = 1500;
	
	private static ArrayList<Short> amplitudes = new ArrayList<Short>();

	public static void main(String[] args) {

		try {

			File myObj = new File("fitness18711.txt");
			Scanner myReader = new Scanner(myObj);
			while (myReader.hasNextLine()) {
				String data = myReader.nextLine();
				if (data.charAt(0) != '(') {
					System.out.println(data);
					amplitudes.add(Short.parseShort(data));
				}
			}
			myReader.close();

			System.out.println("Make sound");
			byte[] buf = new byte[2];

			AudioFormat af = new AudioFormat((float) frequency, 16, 1, true, false);

			SourceDataLine sdl = AudioSystem.getSourceDataLine(af);
			sdl.open();
			sdl.start();


			for (int i = 0; i < amplitudes.size(); i++) { //1000 ms in 1 second
				short a = amplitudes.get(i);

				buf[0] = (byte) (a & 0xFF); //write 8bits ________WWWWWWWW out of 16
				buf[1] = (byte) (a >> 8); //write 8bits WWWWWWWW________ out of 16
				sdl.write(buf, 0, 2);
			}

			sdl.drain();
			sdl.stop();
		} catch (LineUnavailableException e) {
			e.printStackTrace();
		} catch (FileNotFoundException fe) {
			fe.printStackTrace();
		}// try
	}

}
