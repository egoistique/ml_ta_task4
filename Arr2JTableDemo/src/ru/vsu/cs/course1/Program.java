package ru.vsu.cs.course1;

import ru.vsu.cs.util.ArrayUtils;
import ru.vsu.cs.util.SwingUtils;

import java.io.PrintStream;
import java.util.Locale;


public class Program {

    public static class CmdParams {
        public String inputFile;
        public String outputFile;
        public boolean error;
        public boolean help;
        public boolean window;
    }

public static CmdParams parseArgs(String[] args) {
    CmdParams params = new CmdParams();
    int[] pos = {0, 0};

    if (args.length > 0) {
        if (args[0].equals("--help")) {
            params.help = true;
            return params;
        }
        if (args[0].equals("--window")) {
            params.window = true;
            return params;
        }
        if (!args[0].equals("-i") && !args[0].equals("-o")) {
            params.error = true;
            params.help = true;
            return params;
        }
        if (args.length < 2) {
            params.help = true;
            params.error = true;
            return params;
        }

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i")) {
                pos[0] = i + 1;
            }
            if (args[i].equals("-o")) {
                pos[1] = i + 1;
            }
        }
        params.inputFile = args[pos[0]];
        if (args.length > 2) {
            params.outputFile = args[pos[1]];
        }
    } else {
        params.help = true;
        params.error = true;
    }
    return params;
}
    public static void winMain() throws Exception {
        //SwingUtils.setLookAndFeelByName("Windows");
        Locale.setDefault(Locale.ROOT);
        //UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        SwingUtils.setDefaultFont("Microsoft Sans Serif", 18);

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new FrameMain().setVisible(true);
            }
        });
    }

    public static void main(String[] args) throws Exception {
        CmdParams params = parseArgs(args);
        if (params.help) {
            PrintStream out = params.error ? System.err : System.out;
            out.println("Usage:");
            out.println("  <cmd> -i <input-file> (-o <output-file>)");
            out.println("  <cmd> --help");
            out.println("  <cmd> --window  // show window");
            System.exit(params.error ? 1 : 0);
        }
        if (params.window) {
            winMain();
        } else {
            int[][] arr2 = ArrayUtils.readIntArray2FromFile(params.inputFile);
            if (arr2 == null) {
                System.err.printf("Can't read array from \"%s\"%n", params.inputFile);
                System.exit(2);
            }
//            if (params.check) {
//                Task.check(arr2);
//            }

            //Task.check(arr2);
            PrintStream out = (params.outputFile != null) ? new PrintStream(params.outputFile) : System.out;
            out.println(ArrayUtils.toString(arr2));
            out.close();
        }
    }
}
