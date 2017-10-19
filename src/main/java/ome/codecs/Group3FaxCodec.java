/*
 * #%L
 * BSD implementations of Bio-Formats readers and writers
 * %%
 * Copyright (C) 2005 - 2017 Open Microscopy Environment:
 *   - Board of Regents of the University of Wisconsin-Madison
 *   - Glencoe Software, Inc.
 *   - University of Dundee
 * %%
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * #L%
 */

package ome.codecs;

import java.io.EOFException;
import java.io.IOException;
import java.util.Arrays;

import loci.common.DataTools;
import loci.common.RandomAccessInputStream;
import ome.codecs.CodecException;
import ome.codecs.UnsupportedCompressionException;

/**
 */
public class Group3FaxCodec extends BaseCodec {

  // from page 5 of spec
  // the index into the *_TERMINATING_CODES array is the length of the run
  private static final String[] WHITE_TERMINATING_CODES = {
    "00110101",
    "000111",
    "0111",
    "1000",
    "1011",
    "1100",
    "1110",
    "1111",
    "10011",
    "10100",
    "00111",
    "01000",
    "001000",
    "000011",
    "110100",
    "110101",
    "101010",
    "101011",
    "0100111",
    "0001100",
    "0001000",
    "0010111",
    "0000011",
    "0000100",
    "0101000",
    "0101011",
    "0010011",
    "0100100",
    "0011000",
    "00000010",
    "00000011",
    "00011010",
    "00011011",
    "00010010",
    "00010011",
    "00010100",
    "00010101",
    "00010110",
    "00010111",
    "00101000",
    "00101001",
    "00101010",
    "00101011",
    "00101100",
    "00101101",
    "00000100",
    "00000101",
    "00001010",
    "00001011",
    "01010010",
    "01010011",
    "01010100",
    "01010101",
    "00100100",
    "00100101",
    "01011000",
    "01011001",
    "01011010",
    "01011011",
    "01001010",
    "01001011",
    "00110010",
    "00110011",
    "00110100"
  };

  private static final String[] BLACK_TERMINATING_CODES = {
    "0000110111",
    "010",
    "11",
    "10",
    "011",
    "0011",
    "0010",
    "00011",
    "000101",
    "000100",
    "0000100",
    "0000101",
    "0000111",
    "00000100",
    "00000111",
    "000011000",
    "0000010111",
    "0000011000",
    "0000001000",
    "00001100111",
    "00001101000",
    "00001101100",
    "00000110111",
    "00000101000",
    "00000010111",
    "00000011000",
    "000011001010",
    "000011001011",
    "000011001100",
    "000011001101",
    "000001101000",
    "000001101001",
    "000001101010",
    "000001101011",
    "000011010010",
    "000011010011",
    "000011010100",
    "000011010101",
    "000011010110",
    "000011010111",
    "000001101100",
    "000001101101",
    "000011011010",
    "000011011011",
    "000001010100",
    "000001010101",
    "000001010110",
    "000001010111",
    "000001100100",
    "000001100101",
    "000001010010",
    "000001010011",
    "000000100100",
    "000000110111",
    "000000111000",
    "000000100111",
    "000000101000",
    "000001011000",
    "000001011001",
    "000000101011",
    "000000101100",
    "000001011010",
    "000001100110",
    "000001100111"
  };

  // page 7 of the spec
  // index into *_MAKEUP_CODES plus 1 is the multiple of 64
  private static final String[] WHITE_MAKEUP_CODES = {
    "11011",
    "10010",
    "010111",
    "0110111",
    "00110110",
    "00110111",
    "01100100",
    "01100101",
    "01101000",
    "01100111",
    "011001100",
    "011001101",
    "011010010",
    "011010011",
    "011010100",
    "011010101",
    "011010110",
    "011010111",
    "011011000",
    "011011001",
    "011011010",
    "011011011",
    "010011000",
    "010011001",
    "010011010",
    "011000",
    "010011011",
  };

  private static final String[] BLACK_MAKEUP_CODES = {
    "0000001111",
    "000011001000",
    "000011001001",
    "000001011011",
    "000000110011",
    "000000110100",
    "000000110101",
    "0000001101100",
    "0000001101101",
    "0000001001010",
    "0000001001011",
    "0000001001100",
    "0000001001101",
    "0000001110010",
    "0000001110011",
    "0000001110100",
    "0000001110101",
    "0000001110110",
    "0000001110111",
    "0000001010010",
    "0000001010011",
    "0000001010100",
    "0000001010101",
    "0000001011010",
    "0000001011011",
    "0000001100100",
    "0000001100101",
  };

  private static final String EOL = "000000000001";

  // page 8 of spec
  // (index * 64) + 1792 is the run length (color is chosen by subsequent code)
  private static final String[] EXTRA_MAKEUP_CODES = {
    "00000001000",
    "00000001100",
    "00000001101",
    "000000010010",
    "000000010011",
    "000000010100",
    "000000010101",
    "000000010110",
    "000000010111",
    "000000011100",
    "000000011101",
    "000000011110",
    "000000011111",
  };

  @Override
  public byte[] compress(byte[] in, CodecOptions options)
    throws CodecException
  {
    throw new UnsupportedCompressionException("Group 3 Fax Compression not supported");
  }

  @Override
  public byte[] decompress(RandomAccessInputStream in, CodecOptions options)
    throws CodecException, IOException
  {
    if (in == null || in.length() == 0) return null;
    if (options == null) options = CodecOptions.getDefaultOptions();

    byte[] output = new byte[options.maxBytes];

    int pointer = 0;

    boolean isNextRunWhite = true;
    while (pointer < output.length) {
      Run r = getNextCode(in, isNextRunWhite);
      if (r.eol) {
        isNextRunWhite = true;
      }
      else {
        isNextRunWhite = !isNextRunWhite;
      }
      for (int i=0; i<r.length; i++) {
        output[pointer++] = r.value;
      }
    }

    return output;
  }

  private Run getNextCode(RandomAccessInputStream in, boolean checkWhiteTables) throws IOException {
    StringBuffer code = new StringBuffer();

    boolean foundTerminal = false;
    Run run = new Run();
    while (!foundTerminal) {
      int nextBit = in.readBits(1);
      code.append(String.valueOf(nextBit));

      // use endsWith instead of equals to encompass fill bits
      if (code.toString().endsWith(EOL)) {
        run.eol = true;
        break;
      }

      int extraMakeup = DataTools.indexOf(EXTRA_MAKEUP_CODES, code.toString());
      if (extraMakeup >= 0) {
        run.length += ((64 * extraMakeup) + 1792);
        code.setLength(0);
        continue;
      }

      if (checkWhiteTables) {
        int whiteMakeup = DataTools.indexOf(WHITE_MAKEUP_CODES, code.toString());
        if (whiteMakeup >= 0) {
          run.length += (64 * (whiteMakeup + 1));
          code.setLength(0);
          continue;
        }
        int whiteTerminal = DataTools.indexOf(WHITE_TERMINATING_CODES, code.toString());
        if (whiteTerminal >= 0) {
          foundTerminal = true;
          run.length += whiteTerminal;
          run.value = Run.WHITE;
        }
      }
      else {
        int blackMakeup = DataTools.indexOf(BLACK_MAKEUP_CODES, code.toString());
        if (blackMakeup >= 0) {
          run.length += (64 * (blackMakeup + 1));
          code.setLength(0);
          continue;
        }
        int blackTerminal = DataTools.indexOf(BLACK_TERMINATING_CODES, code.toString());
        if (blackTerminal >= 0) {
          foundTerminal = true;
          run.length += blackTerminal;
          run.value = Run.BLACK;
        }
      }
    }
    return run;
  }


  class Run {
    public static final byte WHITE = 0;
    public static final byte BLACK = (byte) 255;

    public int length = 0;
    public byte value = 0x7f;
    public boolean eol = false;

    @Override
    public String toString() {
      return "eol = " + eol + ", length = " + length + ", value = " + value;
    }
  }
}
