using System;
using System.Collections.Generic;
using System.Windows;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Diagnostics;
using System.Text.RegularExpressions;
using System.Threading;
using System.Windows.Controls;
using System.Windows.Forms;
using Button = System.Windows.Controls.Button;
using TextBox = System.Windows.Controls.TextBox;


namespace UAT_EFDC
{

    /// <summary>
    /// MainWindow.xaml contact
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
            
            if (File.Exists("Data.ini"))
            {
                string[] arr1 = File.ReadAllLines("Data.ini");
                if (arr1.Length >= 12)
                {
                    cpt.Text = arr1[0];
                    rft.Text = arr1[1];
                    ift.Text = arr1[2];
                    mpt.Text = arr1[3];
                    gft.Text = arr1[4];
                    mft.Text = arr1[5];
                    t1.Text = arr1[6];
                    t2.Text = arr1[7];
                    t3.Text = arr1[8];
                    t4.Text = arr1[9];
                    t5.Text = arr1[10];
                    ntt.Text = arr1[11];
                }
               
            }
            System.IO.File.Copy(gft.Text, @"getefdc.inp", true);
        }
        private void TextChanged(object sender, TextChangedEventArgs e)
        {
            if (this.IsLoaded)
            {
                string[] arr = new string[12];
                arr[0] = cpt.Text;
                arr[1] = rft.Text;
                arr[2] = ift.Text;
                arr[3] = mpt.Text;
                arr[4] = gft.Text;
                arr[5] = mft.Text;
                arr[6] = t1.Text;
                arr[7] = t2.Text;
                arr[8] = t3.Text;
                arr[9] = t4.Text;
                arr[10] = t5.Text;
                arr[11] = ntt.Text;
                StreamWriter sw = new StreamWriter("Data.ini");
                foreach (string i in arr) { sw.WriteLine(i); }
                sw.Close();
            }
        }
        //import LHS.dll       
        [DllImport("UAT_LIB.dll", CallingConvention = CallingConvention.StdCall)]
        public static extern void LHS();
        private void Sampling()
        {
            
            //1.make file LHSinput.dat=================
            String[] arr = File.ReadAllLines(rft.Text);
            List<string> phead = new List<string> { };
            StreamWriter sw = new StreamWriter(@"LHSinput.dat");
            sw.WriteLine("NOBS  " + t1.Text);
            sw.WriteLine("RANDOM SEED    " + t2.Text);
            foreach (string i in arr)
            {
                if (i != "")
                {
                    string[] arr1 = i.Split(new char[2] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    phead.Add(arr1[0]);
                    sw.WriteLine(arr1[1]);
                    for (int j = 2; j < arr1.Length; j++)
                    {
                        sw.Write(arr1[j] + "       ");
                    }
                    sw.Write("\n");
                }
            }
            sw.Close();
            //2.run LHS sampling===================
            LHS();
            //3.write file Parameters.txt
            StreamWriter sw1 = new StreamWriter(cpt.Text + "\\Parameters.out");
            sw1.WriteLine(string.Join("\t", phead.ToArray()));
            string[] arr2 = File.ReadAllLines("LHSdbg.out");
            foreach (string i in arr2)
            {
                string[] arr3 = i.Split(new char[2] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                List<string> arr4 = new List<string> { };
                for (int j = 2; j < arr3.Length; j++)
                {
                    arr4.Add(Convert.ToDouble(arr3[j]).ToString());
                }
                sw1.WriteLine(string.Join("\t", arr4.ToArray()));
            }
            sw1.Close();
        }
        private void write_inp_Click(object sender, RoutedEventArgs e)
        {
            Sampling();
            //4.Write to .inp
            string[] arr = File.ReadAllLines(cpt.Text + "\\Parameters.out");
            List<string[]> para = new List<string[]> { };
            foreach (string i in arr) { para.Add(i.Split(new char[2] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries)); }
            string[] arr1 = File.ReadAllLines(ift.Text);
            if (Directory.Exists(cpt.Text + "\\inp"))
            {
                Directory.Delete(cpt.Text + "\\inp", true);
            }
            Directory.CreateDirectory(cpt.Text + "\\inp");
            for (int k = 1; k < para.Count; k++)
            {
                for (int i = 1; i < arr1.Length; i++)
                {
                    for (int j = 0; j < para[0].Length; j++)
                    {
                        if (arr1[i].Length >= 2)
                        {
                            if (arr1[i][0] == 'C' && arr1[i][1] != ' ')
                            {
                                int tmp = Array.IndexOf(arr1[i].Split(new char[2] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries), para[0][j]);
                                if (tmp >= 0)
                                {
                                    string[] tmpinp = arr1[i + 1].Split(new char[2] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                                    tmpinp[tmp - 1] = para[k][j];
                                    arr1[i + 1] = "  ";
                                    arr1[i + 1] += string.Join("  ", tmpinp);
                                }
                            }
                        }
                    }
                }
                FileStream fs3 = new FileStream(cpt.Text + "\\inp\\" + k.ToString("0000") + Regex.Match(ift.Text, @"[^/\\]+[/\\]*$"), FileMode.Create, FileAccess.Write);
                StreamWriter sw3 = new StreamWriter(fs3);
                foreach (string i in arr1)
                {
                    sw3.WriteLine(i);
                }
                sw3.Close();
                fs3.Close();
            }
            System.Windows.MessageBox.Show("End of Sampling and Writing!");
        }
        private void run_model_Click(object sender, RoutedEventArgs e)
        {
            if (File.Exists("Settings.ini"))
            {
                string[] settings = File.ReadAllLines("Settings.ini");
                //prepare running model 
               string resstr= cpt.Text + "\\RESULT";
                if(Directory.Exists(resstr))
                { Directory.Delete(resstr, true); }
                    Directory.CreateDirectory(resstr);
                int num = Convert.ToInt32(t3.Text);
                string inpfile = Regex.Match(ift.Text, @"[^/\\]+[/\\]*$").ToString();
                string casepath = cpt.Text;
                System.IO.File.Copy(@"getefdc.inp", mpt.Text + "\\getefdc.inp", true);

                for (int no = 1; no <= num; no++)
                {
                    if (no == 1)
                    {
                        System.IO.File.Copy(mpt.Text + "\\" + inpfile, mpt.Text + "\\" + inpfile + ".backup", true);
                    }
                    System.IO.File.Copy(casepath + "\\inp\\" + no.ToString("0000") + inpfile, mpt.Text + "\\" + inpfile, true);
                    var procInfo = new ProcessStartInfo(settings[0], "-NOP -NT" + ntt.Text)
                    {
                        CreateNoWindow = true,
                        UseShellExecute = true,
                        //WindowStyle = ProcessWindowStyle.Hidden,
                        WorkingDirectory = mpt.Text
                    };
                    Process run = Process.Start(procInfo);
                    run.WaitForExit();
                    
                    var processInfo = new ProcessStartInfo(settings[1], "getefdc.inp")
                    { 
                        CreateNoWindow = true,
                        UseShellExecute = true,
                        //WindowStyle = ProcessWindowStyle.Hidden,
                        WorkingDirectory = mpt.Text
                    };
                    Process outefdc = Process.Start(processInfo);
                    outefdc.WaitForExit();
                    string resstri= cpt.Text + "\\RESULT\\RESULT" + no.ToString("0000");
                    if (Directory.Exists(resstri))
                    { Directory.Delete(resstri, true); }
                    System.IO.Directory.Move(mpt.Text + "\\#output\\RESULT",resstri);
                }
                System.Windows.MessageBox.Show("End of Running!");
            }
            else { System.Windows.MessageBox.Show("No settings!"); }
        }
        //import CALLHV
        [DllImport("UAT_LIB.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int CALLHV(  ref int NO_,  int[] NUM_,   Double[,,] mn_,   Double[,] sc_);
        private void calculate_lhv_Click(object sender, RoutedEventArgs e)
        {
            string sczpath = mft.Text;
            //get scz file
            String[] arr = File.ReadAllLines(sczpath);
            //transfer to arr1
            List<string[]> arr1 = new List<string[]> { };
            foreach (string i in arr) { arr1.Add(i.Split(new char[2] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries)); }
            //get scz matrix
            List<Double[]> scz = new List<Double[]> { };
            foreach (string[] i in arr1)
            {
                if (i[0] != "**")
                {
                    List<Double> tmpd = new List<Double> { };
                    foreach (string j in i) { tmpd.Add(Convert.ToDouble(j)); }
                    scz.Add(tmpd.ToArray()); tmpd.Clear();
                }
            }
            for (int no= Convert.ToInt32(t4.Text); no <= Convert.ToInt32(t5.Text) ; no++)
            {
                //get mnz file
                List<List<Double[]>> mnzs = new List<List<Double[]>> { };
                var files = Directory.GetFiles(cpt.Text+"\\RESULT\\RESULT"+ no.ToString("0000"), svc.SelectionBoxItem.ToString() + "_*_CEL.DAT");
                foreach (string file in files)
                {
                    List<Double[]> mnz = new List<Double[]> { };
                    String[] arrx = File.ReadAllLines(file);
                    //transfer to arr1
                    List<string[]> arrx1 = new List<string[]> { };
                    foreach (string i in arrx) { arrx1.Add(i.Split(new char[2] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries)); }
                    //get mnz matrix
                    foreach (string[] i in arrx1)
                    {
                        foreach (Double[] sci in scz)
                        {
                            if (i[0] != "**" && Convert.ToDouble(i[0]) == sci[0])
                            {
                                List<Double> tmpd = new List<Double> { };
                                foreach (string j in i)
                                {
                                    tmpd.Add(Convert.ToDouble(j));
                                }
                                mnz.Add(tmpd.ToArray());
                                tmpd.Clear();
                                break;
                            }
                        }
                    }
                    mnzs.Add(mnz);
                }
                Double[,,] mn = new Double[mnzs.Count, mnzs[0].Count, mnzs[0][0].Length];
                for (int i = 0; i < mnzs.Count; i++)
                {
                    for (int j = 0; j < mnzs[0].Count; j++)
                    {
                        for (int k = 0; k < mnzs[0][0].Length; k++) { mn[i, j, k] = mnzs[i][j][k]; }
                    }
                }
                Double[,] sc = new Double[scz.Count, scz[0].Length];
                for (int i = 0; i < scz.Count; i++)
                {
                    for (int j = 0; j < scz[0].Length; j++) { sc[i, j] = scz[i][j]; }
                }
                int[] num = new int[] { mn.GetLength(2), mn.GetLength(1), mn.GetLength(0) };
                CALLHV( ref no,  num ,   mn,    sc);
            }
            System.IO.File.Delete(cpt.Text + "\\lhv.out");
            System.IO.File.Move(@"LHV.OUT", cpt.Text + "\\lhv.out");
            System.Windows.MessageBox.Show("Ending of Calculate");
        }
        private void loadsettings(object sender, RoutedEventArgs e)
        {
            Window1 w1 = new Window1();
            w1.Show();
        }

        private void import_dir_Click(object sender, RoutedEventArgs e)
        {
            //get path
            Button b = (Button)sender;
            FolderBrowserDialog dialog = new FolderBrowserDialog();
            dialog.SelectedPath = b.Tag.ToString();
            dialog.Description = "Please select a folder";
            if (dialog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                b.Tag = dialog.SelectedPath;
            }
        }
        private void import_file_Click(object sender, RoutedEventArgs e)
        {
            //get path
            Button b = (Button)sender;
            OpenFileDialog op = new OpenFileDialog();
            op.FileName = b.Tag.ToString();
            op.InitialDirectory=b.Tag.ToString();  
            op.Filter = "inp file(*.inp)|*.inp|all files(*.*)|*.*";
            op.ShowDialog();
            b.Tag = op.FileName;
        }
        private void import_getefdcfile_Click(object sender, RoutedEventArgs e)
        {
            //get path
            import_file_Click(sender,  e);
            System.IO.File.Copy(gft.Text, @"getefdc.inp", true);

        }
        private void loadw2(object sender, RoutedEventArgs e)
        {
           Window2 w2 = new Window2(mpt.Text);
            w2.Show();

        }
    }
}
