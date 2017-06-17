#include <stdio.h>  
#include <omp.h>
/* function: 计算任意位数pi值  
 *    auther: ZhangYachao  
 *       blog:http://blog.csdn.net/u012027907  
 *        */

void CaculatePi()  
{  
	int len,i;                   //len为小数长度  
	int numberator = 1,denominator = 3,result,carry;  
	int flag = 1,count = 0;      //继续循环的标志及循环的次数  
	char *pi,*temp;              //指向保存pi值和临时计算结果的数据  
	printf("请输入小数位数：");  
	scanf("%d",&len);  

	len += 2;   //增加两位  
	if(!(pi = (char*)malloc(sizeof(char)*len)))   //分配保存pi值的内存  
	{  
		printf("分配内存失败！\n");  
		exit(0);  
	}  
	if(!(temp = (char*)malloc(sizeof(char)*len))) //分配保存呢临时值的内存  
	{  
		printf("分配内存失败！\n");  
		exit(0);  
	}  

	for(i = 0; i < len; i++)  //初始化数组  
	{  
		pi[i] = temp[i] = 0;  
	}  
	pi[1] = 2;         //置初值  
	temp[1] = 2;   

	while(flag && (++count < 2147483647))  //int的最大值 2147483647  
	{  
		carry = 0;  
		for(i = len-1; i > 0; i--)     //从低位到高位相乘  
		{  
			result = temp[i] * numberator+carry; //用每一位去乘，再加上进位  
			temp[i] = result % 10;               //保存个数  
			carry = result / 10;                 //进位  
		}  

		carry = 0;  
		for(i = 0; i < len; i++)                 //有高位到低位进行除法运算  
		{  
			result = temp[i] + carry*10;         //当前位加上前一位的余数  
			temp[i] = result / denominator;      //当前位的整数部分  
			carry = result % denominator;        //当前位的余数，累加到下一位的运算  
		}  
		flag = 0; //清除标志  

		for(i = len-1; i > 0; i--)                 
		{     
			result = pi[i] + temp[i];            //将计算结果累加到result中  
			pi[i] = result % 10;                 //保留一位  
			pi[i-1] += result / 10;              //向高位进位  
			flag |= temp[i];                     //若temp中的数全为0,退出循环  
		}   
		numberator++;       //累加分子  
		denominator += 2;   //累加分母  
	}  
	printf("\n计算了%d次\n",count);              //输出循环次数  
	printf("\t--- 第1-1000为小数----\n");  
	printf("PI = \n");  
	printf("%d.",pi[1]);  
	for(i = 2; i < len; i++)  
	{  
		if((i>2) && (i-2)%10 == 0)        //每10位小数间加一个空格  
			printf(" ");  
		if((i>2) && (i-2)%50 == 0)       //每50位小数换行  
			printf("\n");   

		printf("%d",(int)pi[i]);         //输出一位小数  
	}  
	printf("\n");  

}  
void main()  
{ 
	double start_time,end_time;
	start_time=omp_get_wtime();
	CaculatePi();
	end_time=omp_get_wtime();
	printf("all time is %.8fms\n",(end_time-start_time)*1000.0);

} 
