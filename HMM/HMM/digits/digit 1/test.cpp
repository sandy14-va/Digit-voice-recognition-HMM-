#include <iostream>
#include <stdlib.h>

using namespace std;

int main()
{   
    double store[10][5][12],inp[27200],dc_shift=0,max=0.0,min=0.0,temp[320],check=0.0,R[13],C[13],A_f[12];
    double in[7040];
    char c[50];
    int i=0,j=0,k=0;
    FILE * text;
    fopen("digit_1_1.txt","r");
    while(!feof(text))
    {
        fscanf(text,"%s",c);
        int val=atof(c);
        if(val>=1000)
        {
            for(j=0;j<7040;j++)
            {
                fscanf(text,"%s",c);
                in[j]=atof(c)-dc_shift;
            }
            break;
        }
        i++;
    }
    fclose(text);
    cout<<i<<endl;
    for(j=0;j<7040;j++)
    {
        cout<<in[j]<<endl;
    }
}
