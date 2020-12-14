using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Forms;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using Button = System.Windows.Controls.Button;
using TextBox = System.Windows.Controls.TextBox;

namespace UAT_EFDC
{
    /// <summary>
    /// Interaction logic for Window1.xaml
    /// </summary>
    public partial class Window1 : Window
    {
        public Window1()
        {
            InitializeComponent(); 
            if (File.Exists("Settings.ini"))
            {
                string[] arr1 = File.ReadAllLines("Settings.ini");
                if (arr1.Length >= 2)
                {
                    t0.Text = arr1[0];
                    t1.Text = arr1[1];                   
                }
            }
        }

        private void TextChanged(object sender, TextChangedEventArgs e)
        {
            string[] arr = new string[2];
            arr[0] = t0.Text;
            arr[1] = t1.Text;           
            StreamWriter sw = new StreamWriter("Settings.ini");
            foreach (string i in arr) { sw.WriteLine(i); }
            sw.Close();
        }

        private void import_file_Click(object sender, RoutedEventArgs e)
        {
            //get path
            Button b = (Button)sender;
            OpenFileDialog op = new  OpenFileDialog();
            op.FileName = b.Tag.ToString();
            op.FileName = b.Tag.ToString();
            op.Filter = "exe file(*.exe)|*.exe|all files(*.*)|*.*";
            op.ShowDialog();
            b.Tag = op.FileName;
        }
    }
}
