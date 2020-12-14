using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace UAT_EFDC
{
    /// <summary>
    /// Interaction logic for Window2.xaml
    /// </summary>
    public partial class Window2 : Window
    {
        
        public Window2(string mpt)
        {
            InitializeComponent();
            if (File.Exists(@"getefdc.inp"))
            {
               string [] arr = File.ReadAllLines(@"getefdc.inp");
                List<string> al = new List<string> { };
                foreach(string i in arr) { if (i[0] !='*') { al.Add(i); } };
                p.Text = mpt+"\\efdc.inp";              
                al.RemoveAt(0);
               string[] arr1=al[0].Split(new char[2] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);               
                if (arr1.Length >= 9)
                {
                    t1.Text = arr1[0];
                    t2.Text = arr1[1];
                    t3.Text = arr1[2];
                    t4.Text = arr1[3];
                    t5.Text = arr1[4];
                    t6.Text = arr1[5];
                    t7.Text = arr1[6];
                    t8.Text = arr1[7];
                    t9.Text = arr1[8];                   
                }
                al.RemoveAt(0);
                List<string> xyz = new List<string> { };
                foreach (string i in al) { string[] ax = i.Split(new char[2] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries); xyz.Add(string.Join("\t",ax)); }
               
                t10.Text = string.Join("\n",xyz.ToArray());
            }
        }
        private void TextChanged(object sender, TextChangedEventArgs e)
        {
            if (this.IsLoaded)
            {
                string[] arr = new string[3];
                string[] arr1 = new string[9];

                arr[0] = p.Text;
                arr1[0] = t1.Text;
                arr1[1] = t2.Text;
                arr1[2] = t3.Text;
                arr1[3] = t4.Text;
                arr1[4] = t5.Text;
                arr1[5] = t6.Text;
                arr1[6] = t7.Text;
                arr1[7] = t8.Text;
                arr1[8] = t9.Text;
                arr[1] = string.Join("\t", arr1);
                arr[2] = t10.Text;
                StreamWriter sw = new StreamWriter(@"getefdc.inp");
                foreach (string i in arr) { sw.WriteLine(i); }
                sw.Close();
            }
        }

    }
}
