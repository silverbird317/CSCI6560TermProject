import java.util.*;
import java.util.stream.Collectors;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.SourceDataLine;
import javax.sound.sampled.UnsupportedAudioFileException;

public class Body {

	//44100 sample points per 1 second
	private static int frequency = 44100;
	private static int durationMs = 1500;

	private static int populationNum = 20;

	private static Equation[] equations = new Equation[populationNum];
	private static Short[][] population = new Short[populationNum][(int)((float)44100 / 1000 * durationMs)];
	private static int[] fitness = new int[populationNum];

	//private static Equation[] newEQ= new Equation[populationNum];
	//private static Short[][] newGen = new Short[populationNum][(int)((float)44100 / 1000 * durationMs)];
	//private static int[] newFitness = new int[populationNum];

	private static Short[] sampleAmplitudes;
	private static int bucketNum = 57;
	private static int[] sampleAmplitudesBuckets = new int[bucketNum];
	private static int[] sampleBucketDifference = new int[bucketNum - 1];

	private static int maxGeneration = 10000;

	private static int bestFitnessCount = 0;
	private static int bestFitness = Integer.MIN_VALUE;
	private static int same = 0;
	private static int count = 0;
	//private static int fitnessStore;
	//private static Equation equationStore;
	//private static Short[] phenoStore;

	private static short[] dx= new short[(int)((float)44100 / 1000 * durationMs)];

	private static Random random = new Random();

	private static int nodeMax = 50;

	public static void main(String[] args) {
		// get sample recording data
		sampleAmplitudes = getSampleAmplitude("CoughTrimmed.wav");
		sampleAmplitudesBuckets = map(sampleAmplitudes);
		sampleBucketDifference = bucketDifference(sampleAmplitudesBuckets);

		// generate initial population and calculate fitness
		for (int j = 0; j < equations.length; j++) {
			Equation e = new Equation();
			e.generateRandom(nodeMax);
			equations[j] = e;
			//System.out.println("\n" + equations[j].toString());

			for (int i = 0; i < durationMs * (float) 44100 / 1000; i++) { //1000 ms in 1 second

				short a = equations[j].calculate(i);
				population[j][i] = a;
			}

			//Arrays.sort(population[j]);
			int[] phenotype = sortIndividual(population[j], equations[j]);
			int[] indBucketDiff = bucketDifference(phenotype);
			fitness[j] = calculateFitness(indBucketDiff, sampleBucketDifference);
		}
		sortPopulation(0, population.length - 1);

		//Random random = new Random();
		// start evolution
		for (int i = 0; i < maxGeneration; i++) {

			//int[] parentIndices = parentSelection();
			//Equation child = crossover(/*0, 1*/parentIndices[0], parentIndices[1]);

			//System.out.println("Crossovered: " + child + "\nCount: " + child.equationNodeList.size());

			// replace random individual from bottom 20% of population
			//int replace = population.length - random.nextInt((int) (population.length * 0.2)) - 1;
			//equations[replace] = child;
			//population[replace] = generatePhenotype(child);

			//int[] phenotype = sortIndividual(population[replace], equations[replace]);
			//int[] indBucketDiff = bucketDifference(phenotype);
			//fitness[replace] = calculateFitness(indBucketDiff, sampleBucketDifference);

			for (int j = 0; j < population.length; j++) {
				int[] parentIndices = parentSelection();
				Equation child = crossover(/*0, 1*/parentIndices[0], parentIndices[1]);

				child = mutation(child);

				Short[] tempPheno = generatePhenotype(child);
				int[] pheno = sortIndividual(tempPheno, child);
				int[] indBukDiff = bucketDifference(pheno);
				int tempFitness = calculateFitness(indBukDiff, sampleBucketDifference);

				//newEQ[j] = child;
				//newGen[j] = tempPheno;
				//newFitness[j] = tempFitness;

				//newEQ[j] = child;
				//newGen[j] = generatePhenotype(child);
				//int pm = random.nextInt(10);
				//if (pm != 0) {
				//Equation temp = mutation(equations[j]);
				//Short[] tempPheno = generatePhenotype(temp);

				//Equation mutatedKid = mutation(child);
				//System.out.println("    Mutated: " + mutatedKid);

				/*int[] pheno = sortIndividual(tempPheno, temp);
					int[] indBukDiff = bucketDifference(pheno);
					int tempFitness = calculateFitness(indBukDiff, sampleBucketDifference);*/

				if (tempFitness > fitness[j]) {
					equations[j] = child;
					population[j] = tempPheno;
					fitness[j] = tempFitness;
				} // if
				//} // if pm
			} // for

			//equations = newEQ;
			//population = newGen;
			//fitness = newFitness;

			//newEQ= new Equation[populationNum];
			//newGen = new Short[populationNum][(int)((float)44100 / 1000 * durationMs)];
			//newFitness = new int[populationNum];

			sortPopulation(0, population.length - 1);

			if (fitness[0] > bestFitness) {
				bestFitness = fitness[0];
				same = 0;
				writeToFile(population[0], bestFitnessCount);
				bestFitnessCount++;
			} else if (fitness[0] == bestFitness) {
				same++;
			}
			System.out.println(i + ":" + fitness[0]);

			if (same >= 100) {
				same = 0;
				System.out.println("LOCAL OPTIMA, NOAH'S ARK BOTTOM HALF");
				// replace bottom 50% of population and calculate fitness
				for (int j = population.length - 1; j > population.length / 2; j--) {
					Equation e = new Equation();
					e.generateRandom(nodeMax);
					equations[j] = e;
					//System.out.println("\n" + equations[j].toString());

					for (int k = 0; k < durationMs * (float) 44100 / 1000; k++) { //1000 ms in 1 second

						short a = equations[j].calculate(i);
						population[j][k] = a;
					}

					Arrays.sort(population[j]);
					int[] phenotype = sortIndividual(population[j], equations[j]);
					int[] indBucketDiff = bucketDifference(phenotype);
					fitness[j] = calculateFitness(indBucketDiff, sampleBucketDifference);
				}
				sortPopulation(0, population.length - 1);
			} // if
		} // for

		System.out.println("Complete! The best individual was:\n" + equations[0]);
		play(population[0]);

		bestFitnessCount += 10;
		writeToFile(population[1], bestFitnessCount);
		bestFitnessCount++;
		writeToFile(population[2], bestFitnessCount);
		bestFitnessCount++;
		writeToFile(population[3], bestFitnessCount);
		bestFitnessCount++;
		writeToFile(population[4], bestFitnessCount);
	} // main

	private static Short[] generatePhenotype(Equation e) {
		Short[] child = new Short[population[0].length];

		for (int i = 0; i < population[0].length; i++) { //1000 ms in 1 second

			short a = e.calculate(i);
			child[i] = a;
		} // for

		//Arrays.sort(child);

		return child;
	}

	private static void play(Short[] individual) {

		try {
			System.out.println("Make sound");
			byte[] buf = new byte[2];

			AudioFormat af = new AudioFormat((float) frequency, 16, 1, true, false);

			SourceDataLine sdl = AudioSystem.getSourceDataLine(af);
			sdl.open();
			sdl.start();


			for (int i = 0; i < individual.length; i++) { //1000 ms in 1 second
				short a = individual[i];

				buf[0] = (byte) (a & 0xFF); //write 8bits ________WWWWWWWW out of 16
				buf[1] = (byte) (a >> 8); //write 8bits WWWWWWWW________ out of 16
				sdl.write(buf, 0, 2);
			}

			sdl.drain();
			sdl.stop();
		} catch (LineUnavailableException e) {
			e.printStackTrace();
		} // try
	} // play

	/*
	 * Sort population by fitness in decreasing order
	 */
	private static void sortPopulation(int l, int r) {
		if (l < r) {
			// Find the middle point
			int m = l + (r - l) / 2;
			// Sort first and second halves
			sortPopulation(l, m);
			sortPopulation(m + 1, r);
			// Merge the sorted halves
			merge(l, m, r);
		}
	} // sortPopulation

	/*
	 * sortPopulation helper
	 */
	private static void merge(int l, int m, int r) {
		// Find sizes of two subarrays to be merged
		int n1 = m - l + 1;
		int n2 = r - m;

		// Create temp arrays
		int fitnessL[] = new int[n1];
		int fitnessR[] = new int[n2];
		Short popL[][] = new Short[n1][population[0].length];
		Short popR[][] = new Short[n2][population[0].length];

		// Copy data to temp arrays
		for (int i = 0; i < n1; i++) { // left array
			fitnessL[i] = fitness[l + i];
			popL[i] = population[l + i];
		} // for
		for (int j = 0; j < n2; j++) {
			fitnessR[j] = fitness[m + 1 + j];
			popR[j] = population[m + 1 + j];
		} // for

		// Merge the temp arrays

		// Initial indices of first and second subarrays
		int i = 0, j = 0;

		// Initial index of merged subarray array
		int k = l;
		while (i < n1 && j < n2) {
			if (fitnessL[i] >= fitnessR[j]) {
				fitness[k] = fitnessL[i];
				population[k] = popL[i];
				i++;
			} else {
				fitness[k] = fitnessR[j];
				population[k] = popR[j];
				j++;
			} // if-else
			k++;
		} // while

		// Copy remaining elements of L[] if any
		while (i < n1) {
			fitness[k] = fitnessL[i];
			population[k] = popL[i];
			i++;
			k++;
		} // while

		// Copy remaining elements of R[] if any
		while (j < n2) {
			fitness[k] = fitnessR[j];
			population[k] = popR[j];
			j++;
			k++;
		} // while
	} // merge

	private static int[] parentSelection() {
		//Random random = new Random();

		int[] parentPool = new int[8];
		int[] parentFitness = new int[8];
		int[] parents = new int[2];

		// do two tournaments to pick 2 parents
		for (int k = 0; k < 2; k++) {
			// pick 8 potential parents
			for (int i = 0; i < parentPool.length; i++) {
				boolean taken = true;
				int temp = 0;

				while (taken) {
					taken = false;
					temp = random.nextInt(population.length);
					for (int j = 0; j < i; j++) {
						if (parentPool[i] == temp) {
							taken = true;
						} // if
					} // for
				} // while 
				parentPool[i] = temp;
				parentFitness[i] = fitness[parentPool[i]];
			} // for

			// find top parent
			int max = 0;
			for (int i = 1; i < parentPool.length; i++) {
				if (parentFitness[i] > parentFitness[max]) {
					max = i;
				} // if
			} // for - find top parent
			parents[k] = parentPool[max];
		} // for - hold 2 tournaments

		return parents;
	} // parentSelection

	private static Equation crossover(int parent1, int parent2) {
		Equation parentEQ1 = equations[parent1];
		Equation parentEQ2 = equations[parent2];

		Equation child = Equation.crossover(parentEQ1, parentEQ2);

		return child;
	} // crossover

	private static Equation mutation(Equation individual) {

		int mutationType = random.nextInt(8);
		Equation.mutation(individual, mutationType);

		mutationType = random.nextInt(8);
		Equation.mutation(individual, mutationType);
		return individual;
	} // mutation

	private static void survivorSelection() {

	} // survivorSelection

	private static int calculateFitness(int[] indBucketDiff, int[] sampleBucketDiff) {
		int fitness = 0;

		for (int i = 0; i < bucketNum - 1; i++) {
			fitness += Math.abs(indBucketDiff[i] - sampleBucketDiff[i]);
		}

		return -1 * fitness;
	} // calculateFitness

	private static int[] bucketDifference(int[] buckets) {
		int[] bucketDiff = new int[bucketNum - 1];

		for (int i = 0; i < bucketNum - 1; i++) {
			bucketDiff[i] = buckets[i + 1] - buckets[i];
		}

		return bucketDiff;
	} // bucketDifference

	private static void writeToFile(Short[] individual, int index) {
		try {
			HashMap<Short, Integer> sort = new HashMap<Short, Integer>();
			int temp = 0;
			for (int i = 0; i < individual.length; i++) {
				temp = sort.getOrDefault(individual[i], 0);
				if (temp == 0) {
					sort.put(individual[i], 1);
				} else {
					sort.replace(individual[i], temp, temp + 1);
				} // if
			} // for

			Short[] keys = sort.keySet().toArray(new Short[0]);
			Arrays.sort(keys);

			String name = "output" + index + ".txt";

			File myObj = new File(name);
			if (myObj.createNewFile()) {
				System.out.println("File created: " + myObj.getName());
			} else {
				System.out.println("File already exists.");
			}

			FileWriter myWriter = new FileWriter(name);

			myWriter.write(equations[0].toString() + System.lineSeparator());
			for (int i = keys.length - 1; i > - 1; i--) {
				String s = keys[i] + ", " + sort.get(keys[i]);
				myWriter.write(s + System.lineSeparator());
			} // for
			
			String name2 = "amplitudes" + index + ".txt";
			
			File myObj2 = new File(name2);
			if (myObj2.createNewFile()) {
				System.out.println("File created: " + myObj2.getName());
			} else {
				System.out.println("File already exists.");
			}
			
			FileWriter myWriter2 = new FileWriter(name2);

			myWriter2.write(equations[0].toString() + System.lineSeparator());
			for (int i = 0; i < individual.length; i++) {
				String s = individual[i] + "";
				myWriter2.write(s + System.lineSeparator());
			} // for
			
			myWriter.close();
			myWriter2.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static int[] sortIndividual(Short[] dx, Equation eq) {
		try {
			Short[] dy = new Short[dx.length];
			for (int i = 0; i < dy.length; i++) {
				dy[i] = dx[i];
			}
			Arrays.sort(dy);
			int start = 0, end = dy.length;
			//cut off leading and trailing 0s
			for (int i = 0; i < dy.length && dy[i] == 0; i++) {
				start++;
			}
			for (int i = dy.length - 1; i > -1 && dy[i] == 0; i--) {
				end--;
			}

			HashMap<Short, Integer> sort = new HashMap<Short, Integer>();

			//FileWriter myWriter = new FileWriter("output.txt");
			//myWriter.write(eq.toString() + System.lineSeparator());
			/*for (int i = 0; i < dy.length; i++) {
				String s = i + ": " + dy[i];
				myWriter.write(s + System.lineSeparator());
			}*/

			int temp = 0;
			for (int i = 0; i < dy.length; i++) {
				temp = sort.getOrDefault(dy[i], 0);
				if (temp == 0) {
					sort.put(dy[i], 1);
				} else {
					sort.replace(dy[i], temp, temp + 1);
				} // if
			} // for

			Short[] keys = sort.keySet().toArray(new Short[0]);
			Arrays.sort(keys);

			/*for (int i = keys.length - 1; i > - 1; i--) {
				String s = keys[i] + ", " + sort.get(keys[i]);
				myWriter.write(s + System.lineSeparator());
				//System.out.println("key: " + keys[i] + ", value: " + sort.get(keys[i]));
			} // for
			myWriter.close();*/


			// put audio amplitudes in buckets
			int[] store = new int[bucketNum];
			/*
			short min = Short.valueOf(dy[0]);
			short max = Short.valueOf(dy[dy.length - 1]);
			double bucketSize = (max - min) / (double)bucketNum;
			System.out.println("Min: " + min + "\tMax: " + max + "\tBucket Size: " + bucketSize);
			for (int i = 0; i < keys.length; i++) {
				int bucketIndex = ((int)((keys[i] - min) / bucketSize));
				System.out.println(bucketIndex);
				if (bucketIndex == 57) {
					store[56] += sort.get(keys[i]).intValue();
				} else {
					store[bucketIndex] += sort.get(keys[i]).intValue();
				}
			}*/

			store = sortBuckets(dy, keys, sort);

			return store;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	} // sortIndividual

	private static int[] sortBuckets(Short[] arr, Short[] keys, HashMap<Short, Integer> sort) {
		// put audio amplitudes in buckets
		int[] store = new int[bucketNum];

		short min = Short.valueOf(arr[0]);
		short max = Short.valueOf(arr[arr.length - 1]);
		double bucketSize = (max - min) / (double)bucketNum;

		//System.out.println("Min: " + min + "\tMax: " + max + "\tBucket Size: " + bucketSize);
		for (int i = 0; i < keys.length; i++) {
			int bucketIndex = ((int)((keys[i] - min) / bucketSize));
			//System.out.println(bucketIndex);
			if (bucketIndex == 57) {
				store[56] += sort.get(keys[i]).intValue();
			} else {
				store[bucketIndex] += sort.get(keys[i]).intValue();
			}
		} // for
		return store;
	} // sortBuckets

	private static /*<T>*/ int[] map(/*T*/Short arr[]) { // for sample recording
		try {
			int[] store = new int[bucketNum];

			FileWriter myWriter = new FileWriter("CoughTrimmed.txt");

			HashMap</*T*/Short, Integer> sortT = new HashMap</*T*/Short, Integer>();
			/*for (int i = 0; i < arr.length; i++) {
				System.out.println(i + ": " + arr[i]);
			}*/

			int temp = 0;
			for (int i = 0; i < arr.length; i++) {
				//System.out.println(sort.getOrDefault(arr[i], 0));
				temp = sortT.getOrDefault(arr[i], 0);
				if (temp == 0) {
					sortT.put(arr[i], 1);
				} else {
					sortT.replace(arr[i], temp, temp + 1);
				} // if
			} // for

			/*T*/Short[] keys = sortT.keySet().toArray(new Short[0]);
			Arrays.sort(keys);

			for (int i = keys.length - 1; i > - 1; i--) {
				String s = keys[i] + ", " + sortT.get(keys[i]);
				myWriter.write(s + System.lineSeparator());
				//System.out.println("key: " + keys[i] + ", value: " + sortT.get(keys[i]));
			} // for
			myWriter.close();

			/*
			short min = Short.valueOf(arr[0]);
			short max = Short.valueOf(arr[arr.length - 1]);
			double bucketSize = (max - min) / (double)bucketNum;
			System.out.println("Min: " + min + "\tMax: " + max + "\tBucket Size: " + bucketSize);
			for (int i = 0; i < keys.length; i++) {
				int bucketIndex = ((int)((keys[i] - min) / bucketSize));
				//System.out.println(bucketIndex);
				if (bucketIndex == 57) {
					store[56] += sortT.get(keys[i]).intValue();
				} else {
					store[bucketIndex] += sortT.get(keys[i]).intValue();
				}
			}*/

			store = sortBuckets(arr, keys, sortT);
			return store;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	} // map

	private static Short[] getSampleAmplitude(String soundFile) {
		//Wave w1 = new Wave("Cough.wav");
		File file = new File(soundFile);//"CoughTrimmed.wav");
		try {
			AudioInputStream audioInputStream = AudioSystem.getAudioInputStream(file);
			AudioFileFormat format = AudioSystem.getAudioFileFormat(file);
			int numChannels = format.getFormat().getChannels();
			int sampleSizeInBits = format.getFormat().getSampleSizeInBits();
			int frameSize = format.getFormat().getFrameSize();
			float frameRate = format.getFormat().getFrameRate();
			int numBytes = (int) file.length();
			int numFrames = numBytes / frameSize;
			byte[] audioBytes = new byte[numBytes];
			audioInputStream.read(audioBytes);
			Short[] audioData = new Short[numFrames * numChannels];

			for (int i = 0; i < numFrames; i++) {
				for (int j = 0; j < numChannels; j++) {
					int index = (i * numChannels + j) * sampleSizeInBits / 8;
					if (sampleSizeInBits == 16) {
						audioData[i * numChannels + j] = (short) ((audioBytes[index + 1] & 0xff) << 8 | (audioBytes[index] & 0xff));
					} else if (sampleSizeInBits == 8) {
						audioData[i * numChannels + j] = (short) (audioBytes[index] & 0xff);
					} // if
				} // for j
			} // for

			audioInputStream.close();

			int start = 0, end = audioData.length;
			//cut off leading and trailing 0s
			for (int i = 0; i < audioData.length && audioData[i] == 0; i++) {
				start++;
			}
			for (int i = audioData.length - 1; i > -1 && audioData[i] == 0; i--) {
				end--;
			}

			/*Short[] boxedAudioData = new Short[audioData.length];
		    for (int i = 0; i < boxedAudioData.length; i++) {
		    	boxedAudioData[i] = audioBytes[i];
		    }*/

			Short[] trimmedAudioData = Arrays.copyOfRange(audioData, start, end);

			//map(trimmedAudioData);

			//System.out.println(Arrays.toString(audioBytes));
			//System.out.println(Arrays.toString(cutAudioBytes));
			//System.out.println(Arrays.toString(trimmedAudioData));

			Arrays.sort(trimmedAudioData);
			//System.out.println("min: " + trimmedAudioData[0]);
			//System.out.println("max: " + trimmedAudioData[trimmedAudioData.length - 1]);

			return trimmedAudioData;
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		} // try
		return null;		
	}

	static class Equation {
		String equation;
		Node head;
		//int nodeCount;
		//HashMap<Integer, Node> equationNodeList;
		ArrayList<Node> equationNodeList;
		//HashMap<Integer, Integer> sineAmplitudes;
		final static String[] connectingOperations = {"ADD", "SUBTRACT", "MULTIPLY", "DIVIDE"};
		final static String[] subTreeOperations = {"SQUARE", "CUBE"};
		final static String[] sineObject = {"SINE"};
		final static String[] numberObject = {"CONSTANT"};
		Equation(String equation) {
			this.equation = equation;
		} // Equation

		Equation() {
			equationNodeList = new ArrayList<Node>();
			//sineAmplitudes = new HashMap<Integer, Integer>();
			equation = "";
			//nodeCount = 0;
		} // Equation

		public void generateRandom(int nodeMax) {

			//Random random = new Random();
			int randomSinNum = random.nextInt(nodeMax) + 1;

			//System.out.println("max nodes: " + nodeMax + "\tNode Num: " + randomSinNum);

			for (int i = 0; i < randomSinNum; i++) {
				if (i == 0) {
					int tempNumPerSin = random.nextInt(1000);
					Node first = new Node(/*nodeCount, */"SINE", tempNumPerSin);
					head = first;
					equationNodeList.add(first);
					first.key = equationNodeList.indexOf(first);

					// generate random amplitude for sine
					/*int amplitude = random.nextInt(32768);
					sineAmplitudes.put(first.key, amplitude);*/

					i++;
				} else {
					// power or new sine object?
					// 1/5 chance of adding power
					int coinToss = random.nextInt(5);

					if (coinToss == 0) { // power
						int randomSubtree = random.nextInt(equationNodeList.size());

						Node temp = equationNodeList.get(Integer.valueOf(randomSubtree));
						int power = random.nextInt(subTreeOperations.length);
						Node powerNode;
						if (randomSubtree != equationNodeList.indexOf(head)) {
							powerNode = new Node(/*nodeCount, */subTreeOperations[power], temp.parent);

							// replace temp parent's temp child with powerNode
							if (temp.parent.child1 == temp) {
								temp.parent.child1 = powerNode;
							} else {
								temp.parent.child2 = powerNode;
							} // if

							// adopt temp for powerNode and make temp's parent powerNode
							temp.parent = powerNode;
							powerNode.child1 = temp;
						} else {
							powerNode = new Node(/*nodeCount, */subTreeOperations[power]);

							// adopt temp for powerNode and make temp's parent powerNode
							temp.parent = powerNode;
							powerNode.child1 = temp;

							head = powerNode;
						}

						equationNodeList.add(/*powerNode.key, */powerNode);
						powerNode.key = equationNodeList.indexOf(powerNode);
						i++;

					} else { // new sine object and connector
						int randomSubtree = random.nextInt(equationNodeList.size());

						Node part1 = equationNodeList.get(Integer.valueOf(randomSubtree));

						// make other sine node
						int tempNumPerSin = random.nextInt(1000);
						Node part2 = new Node(/*nodeCount, */"SINE", tempNumPerSin);
						i++;
						equationNodeList.add(/*part2.key, */part2);
						part2.key = equationNodeList.indexOf(part2);

						// generate random amplitude for sine
						/*int amplitude = random.nextInt(32768);
						sineAmplitudes.put(part2.key, amplitude);*/

						// choose random connecting operator
						int connectOperator = random.nextInt(connectingOperations.length);
						Node connectNode;
						if (randomSubtree != head.key) {
							connectNode = new Node(/*nodeCount, */connectingOperations[connectOperator], part1.parent);

							// replace sine1 parent's sine1 child with powerNode
							if (part1.parent.child1 == part1) {
								part1.parent.child1 = connectNode;
							} else {
								part1.parent.child2 = connectNode;
							} // if

							// adopt sine1 and sine2 for powerNode and make temp's parent powerNode
							part1.parent = connectNode;
							part2.parent = connectNode;
							connectNode.child1 = part1;
							connectNode.child2 = part2;

						} else {
							connectNode = new Node(/*nodeCount, */connectingOperations[connectOperator]);

							// adopt sine1 and sine2 for powerNode and make temp's parent powerNode
							part1.parent = connectNode;
							part2.parent = connectNode;
							connectNode.child1 = part1;
							connectNode.child2 = part2;

							head = connectNode;
						} // if

						equationNodeList.add(/*connectNode.key, */connectNode);
						connectNode.key = equationNodeList.indexOf(connectNode);
						i++;
					} // if
				} // if
			} // for
		} // generateRandom

		public String toString() {

			String s = toStringHelper(head);
			equation = s;

			return equation;
		} // toString

		private String toStringHelper(Node current) {

			String s = "", s1 = "", s2 = "";

			if (current.child1 != null) {
				s1 = toStringHelper(current.child1);
			} // if

			if (current.child2 != null) {
				s2 = toStringHelper(current.child2);
			} // if

			/*
			 * final static String[] connectingOperations = {"ADD", "SUBTRACT", "MULTIPLY", "DIVIDE"};
		final static String[] subTreeOperations = {"SQUARE", "CUBE"};
		final static String[] sineObject = {"SINE"};
		final static String[] numberObject = {"CONSTANT"};
			 */

			switch (current.nodeType) {
			case "SINE":
				float temp = (float) (current.numberOfSamplesToRepresentFullSin / 2f * Math.PI);
				s += "[ " + current.amplitude + "( SINE( i / " + temp + ") ) ]";
				break;
			case "ADD":
				s += "(" + s1 + " + " + s2;
				break;
			case "SUBTRACT":
				s += "(" + s1 + " - " + s2 + ")";
				break;
			case "MULTIPLY":
				s += "(" + s1 + " * " + s2 + ")";
				break;
			case "DIVIDE":
				s += "(" + s1 + " / " + s2 + ")";
				break;
			case "SQUARE":
				s += "{ " + s1 + " }^2";
				break;
			case "CUBE":
				s += "{ " + s1 + " }^3";
				break;
			} // switch

			//System.out.println(current.key + " ");
			return s;
		} // toStringHelper

		public short calculate(int i) {

			float f = calculateHelper(head, i);

			return (short) f;
		} // toString

		private float calculateHelper(Node current, int i) {

			float f = (float) 0.0, f1 = (float) 0.0, f2 = (float) 0.0;

			if (current.child1 != null) {
				f1 = calculateHelper(current.child1, i);
			} // if

			if (current.child2 != null) {
				f2 = calculateHelper(current.child2, i);
			} // if

			switch (current.nodeType) {
			case "SINE":
				f += current.amplitude * Math.sin(i / (current.numberOfSamplesToRepresentFullSin/ 2.0) * Math.PI);
				//s += "[ " + sineAmplitudes.get(current.key) + "( SINE( i / " + current.numberOfTimesFullSinFuncPerSec + ") ) ]";
				break;
			case "ADD":
				f += f1 + f2;
				//s += s1 + " + " + s2;
				break;
			case "SUBTRACT":
				f += f1 - f2;
				//s += s1 + " - " + s2;
				break;
			case "MULTIPLY":
				f += f1 * f2;
				//s += s1 + " * " + s2;
				break;
			case "DIVIDE":
				f += f1 / f2;
				//s += s1 + " / " + s2;
				break;
			case "SQUARE":
				f += Math.pow(f1, 2);
				//s += "{ " + s1 + " }^2";
				break;
			case "CUBE":
				f += Math.pow(f1, 3);
				//s += "{ " + s1 + " }^3";
				break;
			} // switch

			return f;
		} // calculateHelper

		private Equation deepCopy(Equation parent) {			
			Equation child = new Equation();
			child.head = child.deepCopyHelper(parent.head);

			return child;
		}

		private Node deepCopyHelper(Node current) {
			Node n = null, n1 = null, n2 = null;

			if (current.child1 != null) {
				n1 = deepCopyHelper(current.child1);
			} // if

			if (current.child2 != null) {
				n2 = deepCopyHelper(current.child2);
			} // if

			if (current.nodeType.equals("SINE")) {
				n = new Node(/*nodeCount, */current.nodeType, current.numberOfTimesFullSinFuncPerSec);
				n.amplitude = current.amplitude;
				//sineAmplitudes.put(n.key, amplitude);
			} else {
				n = new Node(/*nodeCount, */current.nodeType);
			} // if

			equationNodeList.add(/*n.key, */n);
			n.key = equationNodeList.indexOf(n);
			//nodeCount++;

			if (n1 != null) {
				n.child1 = n1;
				n1.parent = n;
			}
			if (n2 != null) {
				n.child2 = n2;
				n2.parent = n;
			} // if
			return n;
		} // deepCopyHelper

		private static Equation crossover(Equation parent1, Equation parent2) {

			/*System.out.println("Equation 1: " + parent1);
			System.out.println("Count: " + parent1.nodeCount);
			System.out.println("ArrayList length: " + parent1.equationNodeList.size());
			//System.out.println("Keys: " + parent1.equationNodeList.keySet());
			//parent1.printKeys();
			System.out.println("Equation 2: " + parent2);
			System.out.println("Count: " + parent2.nodeCount);
			System.out.println("ArrayList length: " + parent2.equationNodeList.size());
			//System.out.println("Keys: " + parent2.equationNodeList.keySet());
			//parent2.printKeys();
			 */

			Equation child;

			Equation child1 = new Equation();
			child1 = child1.deepCopy(parent1);

			Equation child2 = new Equation();
			child2 = child2.deepCopy(parent2);

			//System.out.println("Equation 1: " + child1);
			/*System.out.println("Count: " + child1.nodeCount);
			System.out.println("ArrayList length: " + child1.equationNodeList.size());*/
			//System.out.println("Keys: " + child1.equationNodeList.keySet());
			//child1.printKeys();
			//System.out.println("Equation 2: " + child2);
			/*System.out.println("Count: " + child2.nodeCount);
			System.out.println("ArrayList length: " + child1.equationNodeList.size());*/
			//System.out.println("Keys: " + child2.equationNodeList.keySet());
			//child2.printKeys();

			int subtreeIndex1 = random.nextInt(child1.equationNodeList.size());
			Node subtree1 = child1.equationNodeList.get(subtreeIndex1);

			int subtreeIndex2 = random.nextInt(child2.equationNodeList.size());
			Node subtree2 = child2.equationNodeList.get(subtreeIndex2);

			//System.out.println("\nSubtree 1: " + child1.toStringHelper(subtree1));
			//System.out.println("Subtree 2: " + child1.toStringHelper(subtree2));
			/*
			System.out.println("subtree1 " + subtreeIndex1 + " : " + subtree1);
			System.out.println("subtree2 " + subtreeIndex2 + " : " + subtree2);

			System.out.println("Head1? " + (subtree1.key == child1.head.key));
			System.out.println("Head2? " + (subtree2.key == child2.head.key));
			 */

			int child1Length = child1.equationNodeList.size();
			int child2Length = child2.equationNodeList.size();

			//System.out.println("Child1---------------------------------");
			//int subtreeNodeCount1 = countNodesAndAddToList(subtree1, child2Length + 1, child1, child2);
			//System.out.println("Child2---------------------------------");
			//int subtreeNodeCount2 = countNodesAndAddToList(subtree2, child1Length + 1, child2, child1);

			Node temp1 = subtree1.parent;
			Node temp2 = subtree2.parent;
			if (subtreeIndex1 == child1.head.key) {
				//child1.nodeCount = 0;
				child1.equationNodeList = new ArrayList<Node>();
				child1.head = child1.deepCopyHelper(subtree2);
			} else {
				subtree2.parent = temp1;
				if (temp1.child1 == subtree1) {
					temp1.child1 = subtree2;
				} else {
					temp1.child2 = subtree2;
				} // if
			} // if
			if (subtreeIndex2 == child2.head.key) {
				//child2.nodeCount = 0;
				child2.equationNodeList = new ArrayList<Node>();
				child2.head = child2.deepCopyHelper(subtree1);
			} else {
				subtree1.parent = temp2;
				if (temp2.child1 == subtree2) {
					temp2.child1 = subtree1;
				} else {
					temp2.child2 = subtree1;
				} // if
			} // if

			//child1.nodeCount = child1.nodeCount - subtreeNodeCount1 + subtreeNodeCount2;
			reassignKeys(child1);
			child = child1;
			//System.out.println(child1.equationNodeList.size() + " : " + child1.nodeCount);
			child.printKeys();
			//System.out.println("Keys: " + child1.equationNodeList.keySet());
			//child1.printKeys();

			/*System.out.println("Parent 1: " + parent1);
			System.out.println("Parent 2: " + parent2);
			System.out.println(" Child 1: " + child1);
			System.out.println(" Child 2: " + child2);
			 */

			return child;
		} // crossover

		private static int countNodesAndAddToList(Node start, int key, Equation out, Equation in) {
			int count = 0, count1 = 0, count2 = 0;
			if (start.child1 != null) {
				count1 = countNodesAndAddToList(start.child1, key, out, in);
			}
			count += count1;
			key += count1;
			if (start.child2 != null) {
				count2 = countNodesAndAddToList(start.child2, key, out, in);
			}
			count += count2;
			key += count2;
			Node temp = new Node(/*key, */start.nodeType);
			temp.parent = start.parent;
			temp.child1 = start.child1;
			temp.child2 = start.child2;
			temp.amplitude = start.amplitude;
			temp.numberOfSamplesToRepresentFullSin = start.numberOfSamplesToRepresentFullSin;
			temp.numberOfTimesFullSinFuncPerSec = start.numberOfTimesFullSinFuncPerSec;

			System.out.println("add success? " + in.equationNodeList.add(temp));
			temp.key = in.equationNodeList.indexOf(temp);
			System.out.println("remove success? " + out.equationNodeList.remove(start));
			System.out.println("in: " + key);
			System.out.println("out: " + start.key);
			//in.nodeCount++;
			//out.nodeCount--;
			System.out.println("out: " + out.equationNodeList.size());
			System.out.println("in: " + in.equationNodeList.size());
			count += 1;
			return count;
		} // countNodesAndAddToList

		private static void reassignKeys(Equation eq) {
			for (int i = 0; i < eq.equationNodeList.size(); i++) {
				if (i != eq.equationNodeList.get(i).key) {
					eq.equationNodeList.get(i).key = i;
				} // if
			} // for
			/*Set<Integer> keys = eq.equationNodeList.keySet();
			ArrayList<Integer> keysList = new ArrayList<Integer>();
			keysList.addAll(keys);
			Collections.sort(keysList);
			System.out.println("keysList 1: " + keysList);
			ArrayList<Integer> indexList = new ArrayList<Integer>();
			for (int i = 0; i < eq.nodeCount; i++) {
				if (keysList.contains(i)) {
					keysList.remove(keysList.indexOf(i));
				} else {
					indexList.add(i);
				} // if
			} // for
			System.out.println("keysList 2: " + keysList);
			System.out.println("indexList : " + indexList);
			for (int i = 0; i < keysList.size(); i++) {
				Node toRemove = eq.equationNodeList.get(keysList.get(i));
				Node temp = new Node(indexList.get(i), toRemove.nodeType);
				temp.parent = toRemove.parent;
				temp.child1 = toRemove.child1;
				temp.child2 = toRemove.child2;
				temp.amplitude = toRemove.amplitude;
				temp.numberOfSamplesToRepresentFullSin = toRemove.numberOfSamplesToRepresentFullSin;
				temp.numberOfTimesFullSinFuncPerSec = toRemove.numberOfTimesFullSinFuncPerSec;

				eq.equationNodeList.put(indexList.get(i), temp);
				eq.equationNodeList.remove(keysList.get(i));
			} // for
			 */
		} // reassignKeys

		private void removeHelper(Node current) {

			if (current.child1 != null) {
				this.removeHelper(current.child1);
			} // if

			if (current.child2 != null) {
				this.removeHelper(current.child2);
			} // if

			this.equationNodeList.remove(current);
			//this.nodeCount--;
		} // removeHelper

		private void printKeys() {
			equationNodeList.clear();
			//System.out.println("-------------------------");
			printKeyHelper(head);
			//System.out.println("-------------------------");
		}

		private void printKeyHelper(Node current) {
			if (current.child1 != null) {
				printKeyHelper(current.child1);
			}

			equationNodeList.add(current);
			current.key = equationNodeList.indexOf(current);
			//System.out.println("key: " + current.key + ", index: " + equationNodeList.indexOf(current) + ", Node: " + current.nodeType);

			if (current.child2 != null) {
				printKeyHelper(current.child2);
			}
		}

		private static void mutation(Equation parent, int mutationChoice) {
			/*
			 * 1)	Add a random connecting operator and a random sine leaf node to go with it
			 * 2)	Add a random power operator to a random subtree
			 * 3)	Change the connecting operator between two sine expressions
			 * 4)	Change the power operator for a sine expression
			 * 5)	Change the amplitude of a sine leaf node
			 * 6)	Change the angle coefficient for a sine leaf node
			 * 7)	Delete a random connecting operator and its 2 sine leaf children
			 * 8)	Delete a random power operator
			 */
			switch (mutationChoice) {
			case 0:
				addConnectingOperator(parent);
				break;
			case 1:
				addPowerOperator(parent);
				break;
			case 2:
				changeConnectingOperator(parent);
				break;
			case 3:
				changePowerOperator(parent);
				break;
			case 4:
				changeSineAmplitude(parent);
				break;
			case 5:
				changeSineTimesPerSec(parent);
				break;
			case 6:
				deleteConnectingOperator(parent);
				break;
			case 7:
				deletePowerOperator(parent);
				break;
			} // switch
			parent.printKeys();
		} // mutation

		private static void addConnectingOperator(Equation parent) {
			int randomSubtree = random.nextInt(parent.equationNodeList.size());

			Node part1 = parent.equationNodeList.get(Integer.valueOf(randomSubtree));

			// make other sine node
			int tempNumPerSin = random.nextInt(1000) + 1;
			Node part2 = new Node(/*parent.nodeCount, */"SINE", tempNumPerSin);
			//parent.nodeCount++;
			parent.equationNodeList.add(/*part2.key, */part2);
			part2.key = parent.equationNodeList.indexOf(part2);

			// generate random amplitude for sine
			/*int amplitude = random.nextInt(32768);
			sineAmplitudes.put(part2.key, amplitude);*/

			// choose random connecting operator
			int connectOperator = random.nextInt(connectingOperations.length);
			Node connectNode;
			if (randomSubtree != parent.head.key) {
				connectNode = new Node(/*parent.nodeCount, */connectingOperations[connectOperator], part1.parent);

				// replace sine1 parent's sine1 child with powerNode
				if (part1.parent.child1 == part1) {
					part1.parent.child1 = connectNode;
				} else {
					part1.parent.child2 = connectNode;
				} // if

				// adopt sine1 and sine2 for powerNode and make temp's parent powerNode
				part1.parent = connectNode;
				part2.parent = connectNode;
				connectNode.child1 = part1;
				connectNode.child2 = part2;

			} else {
				connectNode = new Node(/*parent.nodeCount, */connectingOperations[connectOperator]);

				// adopt sine1 and sine2 for powerNode and make temp's parent powerNode
				part1.parent = connectNode;
				part2.parent = connectNode;
				connectNode.child1 = part1;
				connectNode.child2 = part2;

				parent.head = connectNode;
			} // if

			parent.equationNodeList.add(/*connectNode.key, */connectNode);
			connectNode.key = parent.equationNodeList.indexOf(connectNode);
			//parent.nodeCount++;
		} // addConnectingOperator

		private static void addPowerOperator(Equation parent) {
			int randomSubtree = random.nextInt(parent.equationNodeList.size());

			Node temp = parent.equationNodeList.get(Integer.valueOf(randomSubtree));
			int power = random.nextInt(subTreeOperations.length);
			Node powerNode;
			if (randomSubtree != parent.head.key) {
				powerNode = new Node(/*parent.nodeCount, */subTreeOperations[power], temp.parent);

				// replace temp parent's temp child with powerNode
				if (temp.parent.child1 == temp) {
					temp.parent.child1 = powerNode;
				} else {
					temp.parent.child2 = powerNode;
				} // if

				// adopt temp for powerNode and make temp's parent powerNode
				temp.parent = powerNode;
				powerNode.child1 = temp;
			} else {
				powerNode = new Node(/*parent.nodeCount, */subTreeOperations[power]);

				// adopt temp for powerNode and make temp's parent powerNode
				temp.parent = powerNode;
				powerNode.child1 = temp;

				parent.head = powerNode;
			}

			parent.equationNodeList.add(/*powerNode.key, */powerNode);
			powerNode.key = parent .equationNodeList.indexOf(powerNode);
			//parent.nodeCount++;
		} // addPowerOperator

		private static void changeConnectingOperator(Equation parent) {
			//System.out.println(parent.nodeCount);
			//System.out.println(parent.equationNodeList.keySet());
			int node = random.nextInt(parent.equationNodeList.size());
			Node temp = parent.equationNodeList.get(node);
			/*System.out.println(node);
			System.out.println(temp);

			System.out.println(parent);

			Integer[] keys = parent.equationNodeList.keySet().toArray(new Integer[0]);
			for (int i = 0; i < keys.length; i++) {
				System.out.print("key: " + keys[i] + "Node: ");
				Node current = parent.equationNodeList.get(keys[i]);
				switch (current.nodeType) {
				case "SINE":
					float temp1 = (float) (current.numberOfSamplesToRepresentFullSin / 2f * Math.PI);
					System.out.println("[ " + current.amplitude + "( SINE( i / " + temp1 + ") ) ]");
					break;
				case "ADD":
					System.out.println(" + ");
					break;
				case "SUBTRACT":
					System.out.println(" - ");
					break;
				case "MULTIPLY":
					System.out.println(" * ");
					break;
				case "DIVIDE":
					System.out.println(" / ");
					break;
				case "SQUARE":
					System.out.println(" ^2 ");
					break;
				case "CUBE":
					System.out.println(" ^3 ");
					break;
				} // switch
			}*/

			for (int i = 0; i < parent.equationNodeList.size() - 1; i++) {
				if (!Arrays.asList(connectingOperations).contains(temp.nodeType)) {
					node = (node + 1) % parent.equationNodeList.size();
					temp = parent.equationNodeList.get(node);
				} else {
					i = parent.equationNodeList.size();
				}
			} // for

			if (!Arrays.asList(connectingOperations).contains(temp.nodeType)) {
				//System.out.println("No connecting operator found");
				return;
			}

			int randomOperator = random.nextInt(connectingOperations.length - 1) + 1;
			int originalIndex = Arrays.asList(connectingOperations).indexOf(temp.nodeType);
			int newIndex = (originalIndex + randomOperator) % connectingOperations.length;
			//System.out.println("original: " + originalIndex);
			//System.out.println("new: " + newIndex);
			temp.nodeType = connectingOperations[newIndex];
		} // changeConnectingOperator

		private static void changePowerOperator(Equation parent) {
			//System.out.println(parent.equationNodeList.size());
			//System.out.println(parent.equationNodeList.keySet());
			int node = random.nextInt(parent.equationNodeList.size());
			Node temp = parent.equationNodeList.get(node);
			/*System.out.println(node);
			System.out.println(temp);

			System.out.println(parent);

			//Integer[] keys = parent.equationNodeList.keySet().toArray(new Integer[0]);
			for (int i = 0; i < parent.equationNodeList.size(); i++) {
				System.out.print("key: " + i + "Node: ");
				Node current = parent.equationNodeList.get(i);
				switch (current.nodeType) {
				case "SINE":
					float temp1 = (float) (current.numberOfSamplesToRepresentFullSin / 2f * Math.PI);
					System.out.println("[ " + current.amplitude + "( SINE( i / " + temp1 + ") ) ]");
					break;
				case "ADD":
					System.out.println(" + ");
					break;
				case "SUBTRACT":
					System.out.println(" - ");
					break;
				case "MULTIPLY":
					System.out.println(" * ");
					break;
				case "DIVIDE":
					System.out.println(" / ");
					break;
				case "SQUARE":
					System.out.println(" ^2 ");
					break;
				case "CUBE":
					System.out.println(" ^3 ");
					break;
				} // switch
			}*/

			for (int i = 0; i < parent.equationNodeList.size() - 1; i++) {
				if (!Arrays.asList(subTreeOperations).contains(temp.nodeType)) {
					node = (node + 1) % parent.equationNodeList.size();
					temp = parent.equationNodeList.get(node);
				} else {
					i = parent.equationNodeList.size();
				}
			} // for

			if (!Arrays.asList(subTreeOperations).contains(temp.nodeType)) {
				//System.out.println("No power operators found");
				return;
			} // if

			int randomOperator = random.nextInt(subTreeOperations.length - 1) + 1;
			int originalIndex = Arrays.asList(subTreeOperations).indexOf(temp.nodeType);
			int newIndex = (originalIndex + randomOperator) % subTreeOperations.length;
			/*System.out.println("key: " + temp.key);
			System.out.println("original: " + originalIndex);
			System.out.println("new: " + newIndex);*/
			temp.nodeType = subTreeOperations[newIndex];
			//System.out.println("new power: " + temp.nodeType);
		} // changePowerOperator

		private static void changeSineAmplitude(Equation parent) {
			int node = random.nextInt(parent.equationNodeList.size());
			Node temp = parent.equationNodeList.get(node);

			/*System.out.println(parent);

			//Integer[] keys = parent.equationNodeList.keySet().toArray(new Integer[0]);
			for (int i = 0; i < parent.equationNodeList.size(); i++) {
				System.out.print("key: " + i + "Node: ");
				Node current = parent.equationNodeList.get(i);
				switch (current.nodeType) {
				case "SINE":
					float temp1 = (float) (current.numberOfSamplesToRepresentFullSin / 2f * Math.PI);
					System.out.println("[ " + current.amplitude + "( SINE( i / " + temp1 + ") ) ]");
					break;
				case "ADD":
					System.out.println(" + ");
					break;
				case "SUBTRACT":
					System.out.println(" - ");
					break;
				case "MULTIPLY":
					System.out.println(" * ");
					break;
				case "DIVIDE":
					System.out.println(" / ");
					break;
				case "SQUARE":
					System.out.println(" ^2 ");
					break;
				case "CUBE":
					System.out.println(" ^3 ");
					break;
				} // switch
			}*/

			while (!temp.nodeType.equals("SINE")) {
				node = (node + 1) % parent.equationNodeList.size();
				temp = parent.equationNodeList.get(node);
			} // while

			//System.out.println("key: " + temp.key);
			temp.amplitude = random.nextInt(32768);
			//System.out.println("new:" + temp.amplitude);
		} // changeSineAmplitude

		private static void changeSineTimesPerSec(Equation parent) {
			int node = random.nextInt(parent.equationNodeList.size());
			Node temp = parent.equationNodeList.get(node);

			/*System.out.println(parent);

			//Integer[] keys = parent.equationNodeList.keySet().toArray(new Integer[0]);
			for (int i = 0; i < parent.equationNodeList.size(); i++) {
				System.out.print("key: " + i + "Node: ");
				Node current = parent.equationNodeList.get(i);
				switch (current.nodeType) {
				case "SINE":
					float temp1 = (float) (current.numberOfSamplesToRepresentFullSin / 2f * Math.PI);
					System.out.println("[ " + current.amplitude + "( SINE( i / " + temp1 + ") ) ]");
					break;
				case "ADD":
					System.out.println(" + ");
					break;
				case "SUBTRACT":
					System.out.println(" - ");
					break;
				case "MULTIPLY":
					System.out.println(" * ");
					break;
				case "DIVIDE":
					System.out.println(" / ");
					break;
				case "SQUARE":
					System.out.println(" ^2 ");
					break;
				case "CUBE":
					System.out.println(" ^3 ");
					break;
				} // switch
			}*/

			while (!temp.nodeType.equals("SINE")) {
				node = (node + 1) % parent.equationNodeList.size();
				temp = parent.equationNodeList.get(node);
			} // while

			//System.out.println("key: " + temp.key);
			temp.numberOfTimesFullSinFuncPerSec = random.nextInt(1000) + 1;
			temp.numberOfSamplesToRepresentFullSin = (float)frequency / temp.numberOfTimesFullSinFuncPerSec;
			//System.out.println("new: " + temp.numberOfSamplesToRepresentFullSin);
		} // changeSineTimesPerSec

		private static void deleteConnectingOperator(Equation parent) {
			// choose a child to delete
			parent.printKeys();
			int sacrifice = random.nextInt(2);

			int node = random.nextInt(parent.equationNodeList.size());
			Node temp = parent.equationNodeList.get(node);

			/*System.out.println(parent);

			//Integer[] keys = parent.equationNodeList.keySet().toArray(new Integer[0]);
			for (int i = 0; i < parent.equationNodeList.size(); i++) {
				System.out.print("key: " + i + "Node: ");
				Node current = parent.equationNodeList.get(i);
				switch (current.nodeType) {
				case "SINE":
					float temp1 = (float) (current.numberOfSamplesToRepresentFullSin / 2f * Math.PI);
					System.out.println("[ " + current.amplitude + "( SINE( i / " + temp1 + ") ) ]");
					break;
				case "ADD":
					System.out.println(" + ");
					break;
				case "SUBTRACT":
					System.out.println(" - ");
					break;
				case "MULTIPLY":
					System.out.println(" * ");
					break;
				case "DIVIDE":
					System.out.println(" / ");
					break;
				case "SQUARE":
					System.out.println(" ^2 ");
					break;
				case "CUBE":
					System.out.println(" ^3 ");
					break;
				} // switch
			}*/

			for (int i = 0; i < parent.equationNodeList.size() - 1; i++) {
				if (!Arrays.asList(connectingOperations).contains(temp.nodeType)) {
					node = (node + 1) % parent.equationNodeList.size();
					temp = parent.equationNodeList.get(node);
				} else {
					i = parent.equationNodeList.size();
				}
			} // for

			if (!Arrays.asList(connectingOperations).contains(temp.nodeType)) {
				//System.out.println("No connecting operators found");
				return;
			} // if

			Node child1 = temp.child1;
			Node child2 = temp.child2;
			Node parentTemp = temp.parent;
			int tempChildNum = -1;
			//System.out.println(temp.key + " vs " + parent.head.key);
			if (temp.key != parent.head.key && parentTemp != null) {
				tempChildNum = (temp == parentTemp.child1) ? 0 : 1;	
			}

			/*System.out.println("key: " + temp.key);
			System.out.println("sacrifice: " + sacrifice);
			System.out.println("child1: " + child1.key);
			System.out.println("child2: " + child2.key);
			if (parentTemp != null) {
				System.out.println("parent: " + parentTemp.key);
			} else {
				System.out.println("parent: " + parentTemp);
			}*/

			Node toRemove;
			Node toKeep;

			if (sacrifice == 0) {
				toRemove = child1;
				toKeep = child2;
			} else {
				toRemove = child2;
				toKeep = child1;
			} // if

			if (tempChildNum == 0) {
				parentTemp.child1 = toKeep;
			} else if (tempChildNum == 1){
				parentTemp.child2 = toKeep;
			} // if
			if (parentTemp == null) {
				parent.head = toKeep;
			}
			toKeep.parent = parentTemp;


			//parent.removeHelper(toRemove);

			//parent.equationNodeList.remove(temp);
			//parent.nodeCount--;
			//System.out.println("keys old: " + parent.equationNodeList.keySet());
			//System.out.println("parent nodeCount: " + parent.equationNodeList.size());
			//reassignKeys(parent);
			//System.out.println("keys: " + parent.equationNodeList.keySet());
		} // deleteConnectingOperator

		private static void deletePowerOperator(Equation parent) {

			/*System.out.println(parent);

			//Integer[] keys = parent.equationNodeList.keySet().toArray(new Integer[0]);
			for (int i = 0; i < parent.equationNodeList.size(); i++) {
				System.out.print("key: " + i + "Node: ");
				Node current = parent.equationNodeList.get(i);
				switch (current.nodeType) {
				case "SINE":
					float temp1 = (float) (current.numberOfSamplesToRepresentFullSin / 2f * Math.PI);
					System.out.println("[ " + current.amplitude + "( SINE( i / " + temp1 + ") ) ]");
					break;
				case "ADD":
					System.out.println(" + ");
					break;
				case "SUBTRACT":
					System.out.println(" - ");
					break;
				case "MULTIPLY":
					System.out.println(" * ");
					break;
				case "DIVIDE":
					System.out.println(" / ");
					break;
				case "SQUARE":
					System.out.println(" ^2 ");
					break;
				case "CUBE":
					System.out.println(" ^3 ");
					break;
				} // switch
			}*/

			int node = random.nextInt(parent.equationNodeList.size());
			Node temp = parent.equationNodeList.get(node);

			for (int i = 0; i < parent.equationNodeList.size() - 1; i++) {
				if (!Arrays.asList(subTreeOperations).contains(temp.nodeType)) {
					node = (node + 1) % parent.equationNodeList.size();
					temp = parent.equationNodeList.get(node);
				} else {
					i = parent.equationNodeList.size();
				}
			} // for

			if (!Arrays.asList(subTreeOperations).contains(temp.nodeType)) {
				//System.out.println("No power operators found");
				return;
			} // if

			Node child1 = temp.child1;
			Node parentTemp = temp.parent;
			int tempChildNum = -1;
			//System.out.println(temp.key + " vs " + parent.head.key);
			if (temp.key != parent.head.key && parentTemp != null) {
				tempChildNum = (temp == parentTemp.child1) ? 0 : 1;	
			} // if

			/*System.out.println("key: " + temp.key);
			System.out.println("child1: " + child1.key);
			if (parentTemp != null) {
				System.out.println("parent: " + parentTemp.key);
			} else {
				System.out.println("parent: " + parentTemp);
			} // if
			 */
			if (parentTemp != null) {
				if (tempChildNum == 0) {
					parentTemp.child1 = child1;
				} else {
					parentTemp.child2 = child1;
				}
			} else {
				parent.head = child1;
			}

			child1.parent = parentTemp;
		} // deletePowerOperator
	} // Equation

	static class Node {
		int key;
		String nodeType;
		Node parent;
		Node child1;
		Node child2;
		float value; // for constants only

		int numberOfTimesFullSinFuncPerSec;
		float numberOfSamplesToRepresentFullSin;
		int amplitude;

		/* 
		 * FOR OPERATORS ONLY 
		 */
		Node (/*int key, */String nodeType) {
			//this.key = key;
			this.nodeType = nodeType;
		} // Node

		Node (/*int key, */String nodeType, Node parent) {
			//this.key = key;
			this.nodeType = nodeType;
			this.parent = parent;
		} // Node

		/*
		 *  FOR SINE ONLY
		 */
		Node (/*int key, */String nodeType, int numberOfTimesFullSinFuncPerSec) {
			//this.key = key;
			this.nodeType = nodeType;
			this.numberOfTimesFullSinFuncPerSec = numberOfTimesFullSinFuncPerSec;
			this.numberOfSamplesToRepresentFullSin = (float)frequency / numberOfTimesFullSinFuncPerSec;
			this.amplitude = random.nextInt(32768);
		} // Node

		Node (/*int key, */String nodeType, int numberOfTimesFullSinFuncPerSec, Node parent) {
			//this.key = key;
			this.nodeType = nodeType;
			this.numberOfSamplesToRepresentFullSin = numberOfTimesFullSinFuncPerSec;
			this.parent = parent;
			this.amplitude = random.nextInt(32768);
		} // Node

		/* 
		 * FOR CONSTANTS ONLY 
		 */
		Node (/*int key, */String nodeType, float value) {
			//this.key = key;
			this.nodeType = nodeType;
			this.value = value;
		} // Node

		Node (/*int key, */String nodeType, Node parent, float value) {
			//this.key = key;
			this.nodeType = nodeType;
			this.parent = parent;
			this.value = value;
		} // Node
	} // Node
} // Body