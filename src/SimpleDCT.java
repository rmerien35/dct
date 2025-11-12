import java.io.*;
import java.lang.Math;
import java.util.*;
import binary.*;

/**
*	Image Compression, DCT / IDCT
*	[YCbCr 4:2:2, DCT 8*8, zigzag]
*
*	@author Ronan Merien <rmerien@hotmail.com>
*
*/
public class SimpleDCT {

	static public int COLS = -1;
	static public int ROWS = -1;
	static public int N = 8;
	static public double quality = 2;

	int[][] imageY;
	int[][] imageCr;
	int[][] imageCb;

	double[][] C = new double[N][N];
	double[][] Ct = new double[N][N];
	double[][] temp = new double[N][N];

	int[][] inputY = new int[N][N];
	int[][] inputCr = new int[N][N];
	int[][] inputCb = new int[N][N];

	int[][] outputY = new int[N][N];
	int[][] outputCr = new int[N][N];
	int[][] outputCb = new int[N][N];


	public void SimpleDCT() {
	}

	protected int round_byte(double a) {
		if (a < 0) return 0;
		else if (a > 255) return 255;
		else return (int) Math.round(a);
	}

	// ----------------------------------------------------------------------------------------------------

	public class Pixel {
		int col, row;

		public Pixel(int i, int j) {
			col = i;
			row = j;
		}
	}

	Pixel[] zigzag8 =
	{
			new Pixel(0, 0),
			new Pixel(0, 1), new Pixel(1, 0),
			new Pixel(2, 0), new Pixel(1, 1), new Pixel(0, 2),
			new Pixel(0, 3), new Pixel(1, 2), new Pixel(2, 1), new Pixel(3, 0),
			new Pixel(4, 0), new Pixel(3, 1), new Pixel(2, 2), new Pixel(1, 3), new Pixel(0, 4), 
			new Pixel(0, 5), new Pixel(1, 4), new Pixel(2, 3), new Pixel(3, 2), new Pixel(4, 1), new Pixel(5, 0), 
			new Pixel(6, 0), new Pixel(5, 1), new Pixel(4, 2), new Pixel(3, 3), new Pixel(2, 4), new Pixel(1, 5), new Pixel(0, 6),
			new Pixel(0, 7), new Pixel(1, 6), new Pixel(2, 5), new Pixel(3, 4), new Pixel(4, 3), new Pixel(5, 2), new Pixel(6, 1), new Pixel(7, 0),
			new Pixel(7, 1), new Pixel(6, 2), new Pixel(5, 3), new Pixel(4, 4), new Pixel(3, 5), new Pixel(2, 6), new Pixel(1, 7), 
			new Pixel(2, 7), new Pixel(3, 6), new Pixel(4, 5), new Pixel(5, 4), new Pixel(6, 3), new Pixel(7, 2),
			new Pixel(7, 3), new Pixel(6, 4), new Pixel(5, 5), new Pixel(4, 6), new Pixel(3, 7),
			new Pixel(4, 7), new Pixel(5, 6), new Pixel(6, 5), new Pixel(7, 4),
			new Pixel(7, 5), new Pixel(6, 6), new Pixel(5, 7),
			new Pixel(6, 7), new Pixel(7, 6),
			new Pixel(7, 7)
	};

	/** 
	* zigzag sequence
	*
	*/

	protected Pixel zigzag(int k) {
		return zigzag8[k];
	}

	// ----------------------------------------------------------------------------------------------------

	double[][] quantumY = 
	{ 
		{ 16, 11, 10, 16, 24, 40, 51, 61 }, 
		{ 12, 12, 14, 19, 26, 58, 60, 55 }, 
		{ 14, 13, 16, 24, 40, 57, 69, 56 }, 
		{ 14, 17, 22, 29, 51, 87, 80, 62 }, 
		{ 18, 22, 37, 56, 68, 109, 103, 77 }, 
		{ 24, 35, 59, 64, 81, 104, 113, 92 }, 
		{ 49, 64, 78, 87, 103, 121, 120, 101 }, 
		{ 72, 92, 95, 98, 112, 100, 103, 99 }
	};

	double[][] quantumCrCb = 
	{ 
		{ 17, 18, 24, 47, 99, 99, 99, 99 }, 
		{ 18, 21, 26, 66, 99, 99, 99, 99 }, 
		{ 24, 26, 56, 99, 99, 99, 99, 99 }, 
		{ 47, 99, 99, 99, 99, 99, 99, 99 }, 
		{ 99, 99, 99, 99, 99, 99, 99, 99 }, 
		{ 99, 99, 99, 99, 99, 99, 99, 99 }, 
		{ 99, 99, 99, 99, 99, 99, 99, 99 }, 
		{ 99, 99, 99, 99, 99, 99, 99, 99 }
	};

	/** 
	* initialyze
	*
	* quantified_DCT[i][j] = DCT[i][j] / quantum[i][j]
	*
	* Cosinus Transform Matrix is C
	* C Transposed Matrix is Ct
	*/

	public void initialyze() 
	{
		int i, j;

		for (j = 0; j < N; j++) {

			C[0][j] = 1.0 / Math.sqrt(N);
			Ct[j][0] = C[0][j];

		}

		for (i = 1; i < N; i++) {
			for (j = 0; j < N; j++) {

				C[i][j] = Math.sqrt(2.0 / N) * Math.cos(((2 * j + 1) * i * Math.PI) / (2.0 * N));
				Ct[j][i] = C[i][j];

			}
		}
	}

	// ----------------------------------------------------------------------------------------------------

	/** 
	* forwardDCT
   *
	* DCT[i][j] = (1/sqrt(2N)) C(x) C(y) SUM(x,0,N-1) SUM(y,0,N-1) pixel[x][y] cos((i*pi*(2x+1))/2N) cos((j*pi*(2y+1))/2N)
	* C(0) = 1/sqrt(2) & C(x) = 1 if x>0
	*
   * DCT = C * image * Ct
	* where temp = image * Ct
	* and   DCT = C * temp
	*	
	* input byte from -128 to +127
	*/

	public void forwardDCT(int[][] image, int[][] DCT) 
	{
		int i, j, k;

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {

				temp[i][j] = 0.0;
				for (k = 0; k < N; k++) temp[i][j] += image[i][k] * Ct[k][j];

			}
		}

		double vtemp;

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)	{

				vtemp = 0.0;
				for (k = 0; k < N; k++)	vtemp += C[i][k] * temp[k][j];

				DCT[i][j] = (int) Math.round(vtemp);

			}
		}
	}

	// ----------------------------------------------------------------------------------------------------

	/** 
	* inverseDCT
	*
	* IDCT[x][y] = (1/sqrt(2N)) SUM(i,0,N-1) SUM(j,0,N-1) C(i) C(j) DCT[i][j] cos((i*pi*(2x+1))/2N) cos((j*pi*(2y+1))/2N)
	* C(0) = 1/sqrt(2) & C(x) = 1 if x>0
	*
	* IDCT = Ct * DCT * C
	* where temp = image * C
	* and IDCT = Ct * temp
	*
	*/

	public void inverseDCT(int[][] DCT, int[][] IDCT)
	{
		int i, j, k;
		double vtemp;

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {

				temp[i][j] = 0.0;
				for (k = 0; k < N; k++) temp[i][j] += DCT[i][k] * C[k][j];

			}
		}

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {

				vtemp = 0.0;
				for (k = 0; k < N; k++)	vtemp += Ct[i][k] * temp[k][j];

				// output byte from -128 to +127
				if (vtemp < -128) 		IDCT[i][j] = -128;
				else if (vtemp > 127) 	IDCT[i][j] = 127;
				else							IDCT[i][j] = (int) Math.round(vtemp);

			}
		}
	}

	

	// ---------------------------------------------------------------------------------
	public void compressFile(String inFile, String outFile) throws Exception {
		try {
			int row, col, i, j, k;
			
			FileInputStream fis = new FileInputStream(inFile);
			BinaryInputStream bis = new BinaryInputStream(new BufferedInputStream(fis));

			FileOutputStream fos = new FileOutputStream(outFile);
			BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));


			// Reading the input BMP imagefile
			int p;

			// BitmapFileHeader (14 bytes)
			
			if ((bis.readByte() != (byte) 'B') || (bis.readByte() != (byte) 'M')) throw new Exception("Not a Bitmap file"); // header = 'BM' (2 bytes)
			p = bis.readBit(32); // BMP file size (4 bytes)
			p = bis.readBit(64); // Reserved & Offset (8 bytes)

			// BitmapInfoHeader (40 bytes)

			p = bis.readBit(32); // info header size = 40 (4 bytes)
			COLS = bis.readBit(32); // width (4 bytes)
			ROWS = bis.readBit(32); // height (4 bytes)
			p = bis.readBit(16); // planes = 1 (2 bytes)
			p = bis.readBit(16); // bitcount = 24 (2 bytes)
			if (p != 24) throw new Exception("Not a 24bits Bitmap file");
			p = bis.readBit(32); // compression = 0 (4 bytes)
			p = bis.readBit(32); // image size (4 bytes)
			p = bis.readBit(32); // parameters (4 bytes)
			p = bis.readBit(32); // parameters (4 bytes)
			p = bis.readBit(32); // parameters (4 bytes)
			p = bis.readBit(32); // parameters (4 bytes)


			bos.writeBit(COLS, 16);
			bos.writeBit(ROWS, 16);

			imageY = new int[ROWS][COLS];
			imageCr = new int[ROWS][COLS];
			imageCb = new int[ROWS][COLS];

			initialyze();

			int blue, green, red, y, cb, cr;

			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {

					blue = bis.readBit(8);
					green = bis.readBit(8);
					red = bis.readBit(8);

					y = (int) (0.299 * red + 0.587 * green + 0.114 * blue);
					cb = (int) (-0.1687 * red - 0.3313 * green + 0.5 * blue); // cb = (int) (0.564*(blue - y));
					cr = (int) (0.5 * red - 0.41874 * green - 0.08130 * blue); // cr = (int) (0.713*(red - y));
					
					imageY[row][col] = y - 128;
					imageCr[row][col] = cr;
					imageCb[row][col] = cb;

				}
			}

			for (row = 0; row < ROWS; row += N) {
				for (col = 0; col < COLS; col += N) {

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {
							inputY[i][j] = imageY[row + i][col + j];
						}
					}

					forwardDCT(inputY, outputY);

	
					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;
						outputY[i][j] = (int) Math.round(outputY[i][j] / quantumY[i][j]);

						bos.writeBit(outputY[i][j] + 1024, 11);
					}
				}
			}

			for (row = 0; row < ROWS; row += 2 * N) {
				for (col = 0; col < COLS; col += 2 * N) {

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {

							inputCr[i][j] = (int) ((imageCr[row	+ 2*i][col + 2*j] 
									+ imageCr[row + 2*i + 1][col + 2*j]
									+ imageCr[row + 2*i][col + 2*j + 1] 
									+ imageCr[row + 2*i + 1][col + 2*j+ 1]) / 4);

							inputCb[i][j] = (int) ((imageCb[row	+ 2*i][col + 2*j] 
									+ imageCb[row + 2*i + 1][col + 2*j]	
									+ imageCb[row + 2*i][col + 2*j + 1]	
									+ imageCb[row + 2*i + 1][col + 2*j + 1]) / 4);

						}
					}

					forwardDCT(inputCr, outputCr);

					forwardDCT(inputCb, outputCb);

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;
						outputCr[i][j] = (int) Math.round(outputCr[i][j] / quantumCrCb[i][j]);

						bos.writeBit(outputCr[i][j] + 1024, 11);
					}

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;
						outputCb[i][j] = (int) Math.round(outputCb[i][j] / quantumCrCb[i][j]);

						bos.writeBit(outputCb[i][j] + 1024, 11);
					}

				}
			}

			bis.close();

			bos.writeEOF();

			bos.flush();
		} catch (EOFException e) {
			System.out.println(e);
		} catch (IOException e) {
			System.out.println("compressFile :" + e);
		}
	}

	// ---------------------------------------------------------------------------------
	public void expandFile(String inFile, String outFile) throws Exception {
		try {
			int row, col, i, j, k;

			FileInputStream fis = new FileInputStream(inFile);
			BinaryInputStream bis =	new BinaryInputStream(new BufferedInputStream(fis));

			FileOutputStream fos = new FileOutputStream(outFile);
			BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));

			byte b = 0;

			COLS = bis.readBit(16);	
			ROWS = bis.readBit(16); 

			imageY = new int[ROWS][COLS];
			imageCr = new int[ROWS][COLS];
			imageCb = new int[ROWS][COLS];

			// 14 bytes
			bos.writeByte((byte) 'B');
			bos.writeByte((byte) 'M');
			bos.writeBit(COLS * ROWS * 3 + 54, 32); // BMP file length

			bos.writeBit(0, 32); // Reserved
			bos.writeByte((byte) 54); // Offset
			bos.writeByte((byte) 0);
			bos.writeByte((byte) 0);
			bos.writeByte((byte) 0);

			// 40 bytes
			bos.writeBit(40, 32); // 40 bytes
			bos.writeBit(COLS, 32); // largeur
			bos.writeBit(ROWS, 32); // hauteur
			bos.writeBit(1, 16);


			bos.writeBit(24, 16); // bits by pixel

			bos.writeBit(0, 32);

			bos.writeBit(COLS * ROWS * 3, 32); // image size

			bos.writeBit(0, 32);
			bos.writeBit(0, 32);
			bos.writeBit(0, 32);
			bos.writeBit(0, 32);


			initialyze();

			int blue, green, red, y, cb, cr;

			for (row = 0; row < ROWS; row += N) {
				for (col = 0; col < COLS; col += N) {


					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;

						inputY[i][j] = bis.readBit(11) - 1024;
						inputY[i][j] = (int) Math.round(inputY[i][j] * quantumY[i][j]);
					}

					inverseDCT(inputY, outputY);

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {
							imageY[row + i][col + j] = outputY[i][j];
						}
					}
				}
			}

			for (row = 0; row < ROWS; row += 2 * N) {
				for (col = 0; col < COLS; col += 2 * N) {

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;

						inputCr[i][j] = bis.readBit(11) - 1024;
						inputCr[i][j] = (int) Math.round(inputCr[i][j] * quantumCrCb[i][j]);
					}

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;

						inputCb[i][j] = bis.readBit(11) - 1024;
						inputCb[i][j] = (int) Math.round(inputCb[i][j] * quantumCrCb[i][j]);
					}

					inverseDCT(inputCr, outputCr);

					inverseDCT(inputCb, outputCb);

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {

							imageCr[row + 2 * i][col + 2 * j] = outputCr[i][j];
							imageCr[row + 2 * i][col + 2 * j + 1] = outputCr[i][j];
							imageCr[row + 2 * i + 1][col + 2 * j] = outputCr[i][j];
							imageCr[row + 2 * i + 1][col + 2 * j + 1] = outputCr[i][j];

							imageCb[row + 2 * i][col + 2 * j] = outputCb[i][j];
							imageCb[row + 2 * i][col + 2 * j + 1] = outputCb[i][j];
							imageCb[row + 2 * i + 1][col + 2 * j] = outputCb[i][j];
							imageCb[row + 2 * i + 1][col + 2 * j + 1] = outputCb[i][j];

						}
					}
				}
			}

			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {

					y = imageY[row][col] + 128;
					cr = imageCr[row][col];
					cb = imageCb[row][col];

					blue = round_byte(y + 1.773 * cb);
					green = round_byte(y - 0.34414 * cb - 0.71414 * cr);
					red = round_byte(y + 1.402 * cr);

					bos.writeBit(blue, 8);
					bos.writeBit(green, 8);
					bos.writeBit(red, 8);
				}
			}

			bos.flush();
			fos.close();
			fis.close();
		} 
		catch (EOFException e) {}
		catch (IOException e) {System.out.println("expandFile :" + e);}
	}

	// ---------------------------------------------------------------------------------------------

	static public void help() 
	{
		System.out.println("SimpleDCT 1.1 powered by Ronan Merien <rmerien@hotmail.com>");
		System.out.println("[YCbCr 4:2:2, DCT 8*8, quantization]");
		System.out.println("Usage: java SimpleDCT [options] imagefile");
		System.out.println("-c		compress bmp_to_dctfile");
		System.out.println("-e		expand dct_to_bmpfile");
		System.out.println();
		//System.out.println("-quality n	quality factor between 1 and 100 (default 20)");
		System.out.println("Examples:  java SimpleDCT -c photo.bmp");
		System.out.println("           java SimpleDCT -e photo.dct");
		System.out.println("");
	}

	static public void main(String args[]) 
	{
		String option = "";
		String imagefile = "";

		int argc = args.length;
		int argIndex = 0;
		String arg;

		while (argc > 0) {
			arg = args[argIndex];

			if (arg.startsWith("-")) {

				if (arg.equals("-quality")) {
					--argc;
					++argIndex;
					arg = args[argIndex];
					int p = Integer.parseInt(arg);
					if ((p >= 1) && (p <= 100)) {
						DCT.quality = p / 10;
					}
				}
				else option = arg;
			
			} 
			else if (imagefile.equals("")) {
				imagefile = arg;
			}

			--argc;
			++argIndex;
		}

		if (imagefile.indexOf(".") != -1) {
			StringTokenizer f = new StringTokenizer(imagefile, ".");
			imagefile = f.nextToken();
		}

		SimpleDCT trans = new SimpleDCT();

		// compress bmp_to_dctfile 
		if (option.equals("-c")) 
		{
			System.out.println("compress " + imagefile + ".bmp to " + imagefile + ".dct");

			try {	
				trans.compressFile(imagefile + ".bmp", imagefile + ".dct"); 
			} 
			catch (Exception e) { System.out.println(e); }

		}

		// expand dct_to_bmpfile
		else if (option.equals("-e")) 
		{
			System.out.println("expand " + imagefile + ".dct to " + imagefile + ".idct.bmp");
			
			try {
				trans.expandFile(imagefile + ".dct", imagefile + ".idct.bmp");
			}
			catch (Exception e) { System.out.println(e); }
			
		}

		// help instructions 
		else {
			help();
		}
	}

}