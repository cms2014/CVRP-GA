/**************************************************************/
/*                   有能力的车辆路径问题(CVRP)                  */
/**************************************************************/
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <string>
#include <cstdio>
using namespace std;

#define READFILENAME "../tc/tai75a.dat"  //数据源文件名
#define WRITEFILENAME "../out/out_tai75a.txt"  //输出文件名

#define popSize 50             // 种群个数
#define MAXGENS 1000           // 最大进化代数
#define PXOVER 0.8             // 交叉概率
#define PMUTATION 0.1         // 变异概率
#define MUTATION_SWAP_NUM 4   // 变异时交换基因的次数
#define PUBLISH_WEIGHT 5000    // 惩罚权重

ifstream in;
ofstream out;

struct cityType
{
	double x;        // 横坐标
	double y;        // 纵坐标
	double demand;   // 需求
};
struct cityType *city;  // 城市
int cityNum;            // 城市数目
double depot_x;         // 仓库横坐标
double depot_y;         // 仓库纵坐标
double capacity;        // 装载能力

int generation;          // 当前代数
int cur_best;            // 最好的个体

struct genoType
{
	int* gene;           // 无重复数字串，表示运输路径
	double fitness;      // 当前个体的适应度
	double rfitness;     // 相对适应度，占总适应度的比例
    double cfitness;     // 累计适应度，用于转轮赌的选择
};

struct genoType *population;     // 当前种群
struct genoType *newPopulation;  // 新种群，用于暂时储存产生的新种群

/** 函数声明 */
void assign(genoType, genoType); 		 // 结构体复制
void generateID(int*, int);  			 // 随机产生个体
void swap(int*, int*);  				 // 交换两个int型的整数
double distance(int*, int, int);         // 计算车辆行驶的距离
double totalCapacity(int *, int , int);  // 计算路线上的总需求
void initialize();    					 // 初始化
void evaluate();      					 // 计算适应值
void keep_the_best();  					 // 保存最好的个体
void select();   						 // 选择
void crossover();  						 // 交叉
void XOver(int, int);  					 // 选择两个父体进行顺序交叉操作
bool isExist(int, int*, int, int); 		 // 判断某值是否存在于另一个数组的某一区段中
void mutate(); 						 	 // 变异
void printID(int *); 					 // 打印个体基因组
void elitist(); 						 // 精英化
void report(); 							 // 输出
void reportRoute(int *); 				 // 根据基因序列输出路线到结果中

/** 函数的实现*/

/**
 * 将genoType类型的b复制给a
 * 参数：a genoType 目标结构体
 *      b genoType 源结构体
 */
void assignGenoType(genoType &a, genoType &b)
{
	a.fitness = b.fitness;
	a.rfitness = b.rfitness;
	a.cfitness = b.cfitness;

	for(int i=0; i<cityNum; ++i)
	{
		a.gene[i] = b.gene[i];
	}
}

/**
 * 随机产生个体，即1~cityNum的无重复数字串
 * 参数：array int* 数组容器
 *      maxNum int 范围上限
 */
void generateID(int *array,int maxNum)
{
	int right;

	for(int i=0; i<maxNum; ++i)
	{
		array[i] = i+1;
	}
	for(int i=0; i<maxNum-1; ++i)
	{
		right = i +1 + rand()%(maxNum-i-1);
		swap(&array[i], &array[right]);
	}
}

/**
 * 交换两个int型的整数
 * 参数：x int* 左交换数
 *      y int* 右交换数
 */
void swap(int *x, int *y)
{
	int temp;

	temp = *x;
	*x = *y;
	*y = temp;
}

/**
 * 计算数组中从 仓库~array[start]~...~array[end]~仓库 的距离，start<=end
 * 参数：array int* 个体的基因组
 *      start int  开始下标
 *      end int    结束下标
 * 返回：double 距离
 */
double distance(int *array, int start, int end)
{
	double dist = 0;
	double pre_x = depot_x;  // 记录上一个点的横坐标
	double pre_y = depot_y;  // 记录上一个点的纵坐标
	double cur_x = 0;  // 当前点的横坐标
	double cur_y = 0;  // 当前点的纵坐标

	for(int i=start; i<=end; ++i)
	{
		cur_x = city[ array[i] ].x;
		cur_y = city[ array[i] ].y;
		dist += sqrt(pow((cur_x- pre_x), 2) + pow((cur_y-pre_y), 2));
		pre_x = cur_x;
		pre_y = cur_y;
	}
	dist += sqrt(pow((depot_x- pre_x), 2) + pow((depot_y-pre_y), 2));
	return dist;
}

/**
 * 计算数组中从 仓库~array[start]~...~array[end]~仓库 的总需求和，start<=end
 * 参数：array int* 个体的基因组
 *      start int  开始下标
 *      end int    结束下标
 * 返回：double 总需求和
 */
double totalCapacity(int *array, int start, int end)
{
	double totalCap = 0;
	for(int i=start; i<=end; ++i)
	{
		totalCap += city[array[i] ].demand;
	}
	return totalCap;
}

/**
 * 读入数据并初始化
 */
void initialize()
{
	double temp;      // 忽略某些输入的一个容器

	in.open(READFILENAME);
	if (!in)
	{
		cout<<">>> 打开文件出错"<<endl;
		exit(1);
	}

	// 读取文件中的数据
	in>>cityNum;
	//cityNum = 10;
	city = new cityType[cityNum+1];
	population = new genoType[popSize+1];
	newPopulation = new genoType[popSize+1];

	in>>temp;        // 忽略目前最优解的读入
	in>>capacity;    // 读入车辆装载能力
	in>>depot_x;     // 读入仓库横坐标
	in>>depot_y;     // 读入仓库纵坐标
	for(int i=1; i<=cityNum; ++i)
	{
		in>>temp;                                   // 忽略城市标号的读入
		in>>city[i].x>>city[i].y>>city[i].demand;   // 读取城市坐标和需求
	}
	in.close();

	// 随机产生初始群体，包括队尾的“最好”的个体，不管会不会产生重复的
	srand((unsigned)time(NULL));
	for(int i=0; i<popSize+1; ++i)
	{
		population[i].fitness = 0;
		population[i].rfitness = 0;
		population[i].cfitness = 0;
		population[i].gene = new int[cityNum];
		generateID(population[i].gene, cityNum);

		newPopulation[i].gene = new int[cityNum];
	}
	// 将案例中最好的个体安插到第一代群体中
	// int bestID[75] = {25,74,75,23,
	// 				20,71,73,72,70,52,53,
	// 				63,64,58,56,62,
	// 				22,1,7,9,10,3,2,11,6,13,
	// 				28,14,51,61,30,45,33,55,46,60,54,
	// 				66,67,26,
	// 				57,39,44,36,
	// 				16,4,5,8,17,
	// 				15,19,21,27,18,24,12,
	// 				65,69,59,68,31,29,42,40,41,34,32,37,48,49,35,43,38,47,50};
	// for(int i=0;i<75;++i)
	// {
	// 	population[0].gene[i] = bestID[i];
	// }

}

/**
 * 计算适应值
 */
void evaluate()
{
	double curCapacity;    // 当前车的负载
	double curCityDemand;  // 当前城市的需求
	double fitness = 0;    // 当前基因的适应值=1/路线长度
	int start;             // 开始运输的城市
	int end;               // 结束运输的城市

	for(int mem=0; mem<popSize; ++mem)
	{
		curCapacity = 0;
		fitness = 0;
		start = 0;
		end = 0;
		for(int i=0; i<cityNum; ++i)
		{
			curCityDemand = city[population[mem].gene[i]].demand;
			curCapacity += curCityDemand;
			if(curCapacity <= capacity)
			{
				end = i;
				if(end == cityNum-1)
				{
					fitness += distance(population[mem].gene, start, end);
				}
			}
			else
			{
				if(start == end) // 如果某条线路只有一个城市，而且车辆还运送不过去就加上惩罚权重，不过在本题的数据里几乎不会出现
				{
					fitness += PUBLISH_WEIGHT;
				}
				fitness += distance(population[mem].gene, start, end);
				start = end = i;
				curCapacity = curCityDemand;
			}
		}
		if(start==cityNum-1)  // 处理边界情况，最后一条线路只有一个城市
		{
			fitness += distance(population[mem].gene, start, end);
			if(curCapacity > capacity)
			{
				fitness += PUBLISH_WEIGHT;
			}
		}
		population[mem].fitness = 1.0/fitness;
	}
}

/**
 * 保存最好的个体
 */
void keep_the_best()
{
	cur_best = 0;
	for(int mem=0; mem<popSize; ++mem)
	{
		if (population[mem].fitness > population[popSize].fitness)
        {
            cur_best = mem;
            population[popSize].fitness = population[mem].fitness;
      	}
	}
	for(int i=0; i<cityNum; ++i)
	{
		population[popSize].gene[i] = population[cur_best].gene[i];
	}
}

/**
 * 选择
 */
void select()
{
	int mem, i, j;
	double sum = 0;
	double p;

	// 计算适应值之和
	for(mem=0; mem<popSize; ++mem)
	{
		sum += population[mem].fitness;
	}

	// 计算相对适应值
	for(mem=0; mem<popSize; ++mem)
	{
		population[mem].rfitness = population[mem].fitness / sum;
	}

	// 计算累积适应值，用于转盘赌选择
	population[0].cfitness = population[0].rfitness;
	for(mem=1; mem<popSize; ++mem)
	{
		population[mem].cfitness = population[mem-1].cfitness + population[mem].rfitness;
	}

	// 转盘赌选择
	for(i=0; i<popSize; ++i)
	{
		p = rand()%1000/1000.0;
		if(p < population[0].cfitness)  // 如果概率落在第一个个体
			assignGenoType(newPopulation[i], population[0]);
		else
		{
			for(j=0; j<popSize-1; ++j)
			{
				if(p >= population[j].cfitness && p < population[j+1].cfitness)
					assignGenoType(newPopulation[i], population[j+1]);
			}
		}
	}
	for(i=0; i<popSize; ++i)
	{
		assignGenoType(population[i], newPopulation[i]);
	}

}

/**
 * 交叉
 */
void crossover()
{
	int one;  //选中的；
	int first = 0; // 计数器，用于选择两个父体
	double x;

	for(int mem=0; mem<popSize; ++mem)
	{
		x = rand()%1000/1000.0;
		if(x < PXOVER)
		{
			++first;
			if(first % 2 == 0)
			{
				XOver(one, mem);
			}
			else
				one = mem;
		}
	}
}

/**
 * 选择两个父类进行交叉操作(顺序交叉法)
 * 参数：one int 第一个父类的index
 *      two int 第二个父类的index
 */
void XOver(int one, int two)
{
	// cout<<"父类1="<<endl;
	// printID(population[one].gene);
	// cout<<"父类2="<<endl;
	// printID(population[two].gene);

	int i,j;
	int *newOne = new int[cityNum];
	int *newTwo = new int[cityNum];
	int start = rand() % cityNum;
	int end = rand() % cityNum;
	if(start > end)
		start = end;

	// cout<<"start="<<start<<"，end="<<end<<endl;

	memset(newOne, 0, cityNum);
	memset(newTwo, 0, cityNum);
	for(i=start; i<=end; i++)
	{
		newOne[i] = population[one].gene[i];
		newTwo[i] = population[two].gene[i];
	}

	// 从第二个父个体顺序选择剩下的填充到第一个子类的个体中
	for(i=0, j=0; i<cityNum && j<cityNum; )
	{
		if(i >= start && i <= end)
		{
			++i;
			continue;
		}
		if(!isExist(population[two].gene[j], newOne, start, end))
		{
			newOne[i] = population[two].gene[j];
			++i;
		}
		++j;
	}

	// 从第一个父个体顺序选择剩下的填充到第二个子类的个体中
	for(i=0, j=0; i<cityNum && j<cityNum; )
	{
		if(i >= start && i <= end)
		{
			++i;
			continue;
		}
		if(!isExist(population[one].gene[j], newTwo, start, end))
		{
			newTwo[i] = population[one].gene[j];
			++i;
		}
		++j;
	}

	// 将杂交后产生的新个体复制回父个体中
	for(i=0; i<cityNum; ++i)
	{
		population[one].gene[i] = newOne[i];
		population[two].gene[i] = newTwo[i];
	}

	delete [] newOne;
	delete [] newTwo;
}

/**
 * 判断val是否存在从array[start]到array[end]的数种
 * 参数：val int 要查询的值
 *      array int* 数组
 *      start int 起点下标
 *      end int 终点下标
 * 返回：bool 是否存在
 */
bool isExist(int val, int *array, int start, int end)
{
	for(int i=start; i<=end; ++i)
	{
		if(val == array[i])
			return true;
	}
	return false;
}

/**
 * 变异
 */
void mutate()
{
	int mem, one, two;
	double p;
	for(mem=0; mem<popSize; ++mem)
	{
		p = rand()%1000/1000.0;
		if(p < PMUTATION)
		{
			for(int count=0; count<MUTATION_SWAP_NUM; ++count)
			{
				// do-while 循环保证两个基因不一样
				do {
					one = rand()%cityNum;
					two = rand()%cityNum;
				}while(one==two && cityNum>1);
				swap(&population[mem].gene[one], &population[mem].gene[two]);
			}
		}
	}

}

/**
 * 精英化，popilation[popSize]已经储存了上一代最好的个体
 */
void elitist()
{
	int i;
	int best_mem = 0, worst_mem = 0;
	double best = population[0].fitness;
	double worst = population[0].fitness;

	// 选出当前群里中最好的和最差的个体
	for(i=1; i<popSize; ++i)
	{
		if(population[i].fitness > best)
		{
			best = population[i].fitness;
			best_mem = i;
		}
		else if(population[i].fitness < worst)
		{
			worst = population[i].fitness;
			worst_mem = i;
		}
	}

	// 如果当前最好的个体比上一代最好的个体还好，那么就保存到数组的最后
	// 否则将当前最差的个体替换成上一代最好的个体
	if(best > population[popSize].fitness)
		assignGenoType(population[popSize], population[best_mem]);
	else
		assignGenoType(population[worst_mem], population[popSize]);

}


void printID(int *gene)
{
	for(int i=0; i<cityNum; ++i)
	{
		out<<gene[i]<<" ";
	}
	out<<endl;
}

void report()
{
	out<<"第【"<<generation<<"】代，最好的个体总行驶距离为【"<<1.0/population[popSize].fitness<<"】，基因组序列为："<<endl;
	printID(population[popSize].gene);
}

void reportRoute(int *gene)
{
	double curCapacity = 0;    // 当前车的负载
	double curCityDemand = 0;  // 当前城市的需求
	int start = 0;             // 开始运输的城市
	int end = 0;               // 结束运输的城市
	int i, j;

	out<<"该个体所对应的路线如下："<<endl;
	for(i=0; i<cityNum; ++i)
	{
		curCityDemand = city[gene[i]].demand;
		curCapacity += curCityDemand;
		if(curCapacity <= capacity)
		{
			end = i;
			if(end == cityNum-1)
			{
				for(j=start; j<=end; ++j)
				{
					out<<gene[j]<<" ";
				}
				out<<"| "<<distance(gene, start, end)<<" "<<curCapacity<<endl;
			}
		}
		else
		{
			for(j=start; j<=end; ++j)
			{
				out<<gene[j]<<" ";
			}
			out<<"| "<<distance(gene, start, end)<<" "<<curCapacity-curCityDemand<<endl;
			start = end = i;
			curCapacity = curCityDemand;
		}
	}
	if(start==cityNum-1)  // 处理边界情况，最后一条线路只有一个城市
	{
		for(j=start; j<=end; ++j)
		{
			out<<gene[j]<<" ";
		}
		out<<"| "<<distance(gene, start, end)<<" "<<curCapacity<<endl;
	}
}

int main()
{
	clock_t t_start, t_finish;
	t_start = clock();
	cout<<">>> 《基于遗传算法解决CVRP问题》"<<endl;
	cout<<">>> 遗传算法开始"<<endl;
	out.open(WRITEFILENAME);
	if(!out)
	{
		cout<<">>> 写入文件打开出错！"<<endl;
		return 0;
	}
	generation = 0;
	cout<<">>> 读取数据中......"<<endl;
	initialize();
	cout<<">>> 数据读取成功，初代群体初始化完成"<<endl;
	evaluate();
	keep_the_best();
	report();
	cout<<">>> 群体进化中......"<<endl;
	while(generation < MAXGENS)
	{
		generation++;
		select();
      	crossover();
      	mutate();
      	evaluate();
      	elitist();
      	report();
	}
	out<<endl<<"【结果】"<<endl;
	out<<"经过【"<<MAXGENS<<"】代进化，最好的个体总行驶距离为【"<<1.0/population[popSize].fitness<<"】，基因组序列为："<<endl;
	printID(population[popSize].gene);
	reportRoute(population[popSize].gene);
	out.close();

	t_finish = clock();
	cout<<">>> 群体进化终止，请到out目录下查看过程日志和结果"<<endl
		<<">>> 遗传算法完成，用时"<<(double)(t_finish-t_start)/CLOCKS_PER_SEC<<"秒"<<endl
		<<">>> 版权所有，违版必究"<<endl;
	return 0;
}
