#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include<sys/time.h>
#include <stdint.h>
#include <pthread.h>
#include <fcntl.h>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <unordered_map>
using namespace std;

#define LL long long
#define US unsigned short
#define Max_ID 10000000 //最大结点ID小于这个数用数组做映射
#define MAX_LOG 2000005 // 输入最大边数
#define MAX_circle 20000000 //输出最大环数
#define NT_F 4  // 找环线程数
#define block_size 256 // 找环负载均衡时每一片大小
#define path_size 200 // 反向找三层路径时每个结点的最大路径数
const uint8_t T = 30; // 划分菊花点/非菊花点超参数

const char* input_path = "/data/test_data.txt";
const char* output_path = "/projects/student/result.txt";

struct Timer
{
    /*
        计时函数，可以用于多线程内精准计时
    */
    timeval tic, toc;

    Timer()
    {
        gettimeofday(&tic,NULL);
    }

    void stop(const char* name)
    {
        gettimeofday(&toc,NULL);
        printf("%s: %f(s)\n", name, float(toc.tv_sec-tic.tv_sec) + float(toc.tv_usec-tic.tv_usec)/1000000);
    }
};

struct Edge
{
    int v, w;
};

struct Edge_
{
    int u,v,w;
};
int l1[2*MAX_LOG], l1_len;

struct Node
{
    int n1, n2, wtail, whead;
};

// 读入数据块
struct DBlock
{
    int begin, end, log_len;
    Edge_ *Log;
    int *in_du, *out_du;
};

// 排序数据块
struct SBlock
{
    int begin, end;
};

// 找环数据块
struct CBlock
{
    int n_id;                                 // 分配的ID数
    int *ids;                                 // 分配的ID
    int ns[5];                                // 每种环的个数
    int *ccs[5];                              // 结果
    int *ncs[5];                              // 每个ID每种环的个数
    int *ics[5];                              // 每个ID每种环在cc中的位置
};

// 转换存储范围
struct CRange  // i, [begin, end)
{
    int i, begin, end;
};
int ifsize;                         // 输入文件大小
int max_id, max_edge;
int* ids;                           // 原始ID
Edge *G, *Gb;                       // 非菊花图图，正向反向邻接表
Edge *denseG, *denseGb;            // 菊花图，正向反向邻接表
int *gidx, *gbidx;                  // ID在图中的位置索引
uint8_t *ideg, *odeg;                // 出入度
US *denseIdeg, *denseOdeg;          // 菊花图出入度
bool *isDenseG, *isDenseGb;         // 是否是正向，反向图的菊花点

unordered_map<int, int>id2numMap;
int id2num[Max_ID];
char* node2str;                    // 映射字符串
char* node2len;                     // 映射字符串长度

SBlock sblock[4];                   // 排序数据块
DBlock dblock[4];                   // 读入数据块
CBlock cblock[NT_F];                // 找环数据块
CRange crange[8];                   // 转换存储范围
volatile bool done[8];           // 线程完成标识
volatile int fsize;                 // 输出文件大小
int fd;                             // 输出文件句柄
char *ifmap;                        // mmap

volatile int task;                  // 当前任务ID
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER; //负载均衡用mtx
pthread_cond_t  cond = PTHREAD_COND_INITIALIZER;

int each_r[5] = {0};                // 每种环数量
int total_r = 0;                    // 环路总数
int **pos;                          // 所有环路的指针
int *plen;                          // 所有环路对应长度

int get_file_size(const char *path)
{
    /*
        获取文件大小
    */
    struct stat statbuff;
    stat(path, &statbuff);
    return statbuff.st_size;
}

// 读数据
void *read(void* arg)
{
    /*
        多线程mmap读入
    */
    DBlock* blk = (DBlock*)arg;
    Edge_ *Log = blk->Log;
    char c;
    int i = blk->begin, j = 0;
    int node1, node2, node3;
    while (i < blk->end)
    {
        node1 = 0, node2 = 0, node3 = 0;
        while ((c = ifmap[i++]) != ',')  node1 = (node1 * 10 + c - 48);
        while ((c = ifmap[i++]) != ',')  node2 = (node2 * 10 + c - 48);
        while ((c = ifmap[i++]) >= 48) node3 = (node3 * 10 + c - 48);
        if (ifmap[i] < 48) i++;
        Log[j].u = node1, Log[j].v = node2, Log[j++].w = node3;
    }
    blk->log_len = j;
}

void *arrayHash(void* arg)
{
    /*
        数组映射
    */
    DBlock* blk = (DBlock*)arg;
    Edge_* Log = blk->Log;
    blk->in_du = (int*)calloc(max_id, sizeof(int));
    blk->out_du = (int*)calloc(max_id, sizeof(int));
    int *in_du = blk->in_du, *out_du = blk->out_du;
    int len = blk->log_len;
    for (int j = 0; j < len; j++)
    {
        Edge_ temp = Log[j];
        int u = id2num[temp.u], v = id2num[temp.v];
        in_du[v]++, out_du[u]++;
        Log[j].u = u, Log[j].v = v;
    }
}

void *mapHash(void* arg)
{
    /*
        map映射
    */
    DBlock* blk = (DBlock*)arg;
    Edge_* Log = blk->Log;
    blk->in_du = (int*)calloc(max_id, sizeof(int));
    blk->out_du = (int*)calloc(max_id, sizeof(int));
    int *in_du = blk->in_du, *out_du = blk->out_du;
    int len = blk->log_len;
    for (int j = 0; j < len; j++)
    {
        Edge_ temp = Log[j];
        int u = id2numMap[temp.u], v = id2numMap[temp.v];
        in_du[v]++, out_du[u]++;
        Log[j].u = u, Log[j].v = v;
    }
}

void itoa(int d, char* c, char* l)
{
    /*
        将结点转为字符串
    */
    char* s = (char*)malloc(11*sizeof(char));
    int i = 10;
    s[i--] = ',';
    while (d > 0)
    {
        s[i--] = d % 10 + 48;
        d /= 10;
    }
    *l = 10 - i;
    memcpy(c, s+i+1, 10 - i);
}

void *get_strmap(void* arg)
{
    /*
        将结点转为字符串
        为了保证地址尽量连续，提高cache命中率，node2str数组使用一维指针
    */
    node2str = (char*)malloc(max_id * 16 * sizeof(char));
    node2len = (char*)malloc(max_id * sizeof(char));
    if (l1[0] == 0){
        node2str[0] = '0'; node2str[1] = ','; node2len[0] = 2;
        for (int i = 1; i < max_id; i++)  itoa(l1[i], node2str + i * 16, node2len + i);
    }
    else{
        for (int i = 0; i < max_id; i++)  itoa(l1[i], node2str + i * 16, node2len + i);
    }
}

bool cmp(Edge a, Edge b)
{
    return a.v < b.v;
}

void *sort_edge(void* arg)
{
    /*
        排序
    */
    SBlock* blk = (SBlock*)arg;
    Edge *p1, *p2;
    int begin = blk->begin, end = blk->end;
    for (int i = begin; i < end; i++)
    {
        if (isDenseG[i]){
            p1 = denseG + gidx[i];
            sort(p1, p1 + denseOdeg[i], cmp);
        }
        else{
            p1 = G + i * T;
            sort(p1, p1 + odeg[i], cmp);
        }
    }
    for (int i = begin; i < end; i++)
    {
        if (isDenseGb[i]){
            p2 = denseGb + gbidx[i];
            sort(p2, p2 + denseIdeg[i], cmp);
        }
        else{
            p2 = Gb + i * T;
            sort(p2, p2 + ideg[i], cmp);
        }
    }
}

void* create_graph(void* arg)
{
    /*
        建图
        针对菊花点和非菊花点构建不同的图
        多线程对正、反四个图排序
    */
    // 初始化
    G = (Edge*)malloc(max_id * T * sizeof(Edge));
    Gb = (Edge*)malloc(max_id * T * sizeof(Edge));
    ideg = (uint8_t*)calloc(max_id, sizeof(uint8_t));
    odeg = (uint8_t*)calloc(max_id, sizeof(uint8_t));
    isDenseGb = (bool*)calloc(max_id, sizeof(bool));
    isDenseG = (bool*)calloc(max_id, sizeof(bool));
    denseIdeg = (US*)calloc(max_id, sizeof(US));
    denseOdeg = (US*)calloc(max_id, sizeof(US));
    US *tmp_ideg = (US*)calloc(max_id, sizeof(US));
    US *tmp_odeg = (US*)calloc(max_id, sizeof(US));
    for(int i = 0; i < 4; i++){
        int *in = (dblock + i)->in_du, *out = (dblock + i)->out_du;
        for(int j = 0; j < max_id; j++) tmp_ideg[j] += in[j], tmp_odeg[j] += out[j];
    }
    // 分别菊花点和非菊花点
    gidx = (int*)malloc((max_id + 1) * sizeof(int));
    gbidx = (int*)malloc((max_id + 1) * sizeof(int));
    gidx[0] = 0, gbidx[0] = 0;
    for(int i = 0; i < max_id; i++){
        if (tmp_odeg[i] > T) gidx[i + 1] = gidx[i] + tmp_odeg[i], isDenseG[i] = 1, odeg[i] = tmp_odeg[i], denseOdeg[i] = tmp_odeg[i];
        else gidx[i + 1] = gidx[i], odeg[i] = tmp_odeg[i];
        if (tmp_ideg[i] > T) gbidx[i + 1] = gbidx[i] + tmp_ideg[i], isDenseGb[i] = 1, denseIdeg[i] = tmp_ideg[i];
        else gbidx[i + 1] = gbidx[i], ideg[i] = tmp_ideg[i];
    }
    denseG = (Edge*)malloc(gidx[max_id] * sizeof(Edge));
    denseGb = (Edge*)malloc(gbidx[max_id] * sizeof(Edge));

    // 3. 构建两种邻接表
    US *tmp = (US*)calloc(max_id, sizeof(US));
    US *tmpb = (US*)calloc(max_id, sizeof(US));
    int p1, p2;
    for (int i = 0; i < 4; i++)
    {
        Edge_ *Log = dblock[i].Log;
        int len = dblock[i].log_len;
        for(int j = 0; j < len; j++)
        {
            int u = Log[j].u, v = Log[j].v, w = Log[j].w;
            Edge e1 = {u, w}, e2 = {v, w};
            if (isDenseG[u]){
                p1 = gidx[u];
                denseG[p1 + (tmp[u]++)] = e2;
            }
            else{
                p1 = u * T;
                G[p1 + (tmp[u]++)] = e2;
            }
            if (isDenseGb[v]){
                p2 = gbidx[v];
                denseGb[p2 + (tmpb[v]++)] = e1;
            }
            else{
                p2 = v * T;
                Gb[p2 + (tmpb[v]++)] = e1;
            }
        }
    }
    sblock[0].begin = 0;
    for(int i = 1; i < 4; i++) sblock[i].begin = i * max_id / 4;
    for(int i = 0; i < 3; i++) sblock[i].end = sblock[i+1].begin;
    sblock[3].end = max_id;
    pthread_t thread[4];
    for (int i = 0; i < 4; i++)  pthread_create(&thread[i], NULL, sort_edge, sblock + i);
    for (int i = 0; i < 4; i++)  pthread_join(thread[i], NULL);
}

void load()
{
    /*
        读入文件并建图
        多线程mmap读入文件
        判断最大结点大小，若小于1e7则使用数组做hash，否则使用unordered_map
    */
    Timer timer;
    // 确定每个线程读文件区间
    ifsize = get_file_size(input_path);
    int ifd = open(input_path, O_RDONLY, S_IRUSR);
    ifmap = (char *)mmap(NULL, ifsize, PROT_READ, MAP_SHARED, ifd, 0);
    int dsz = ifsize / 4;
    dblock[0].begin = 0;
    for (int i = 1; i < 4; i++)
    {
        int rp = dsz * i;
        while (ifmap[rp++] != '\n');
        dblock[i].begin = rp;
    }
    for (int i = 0; i < 3; i++)  dblock[i].end = dblock[i+1].begin;
    dblock[3].end = ifsize;
    // 初始化变量
    for (int i = 0; i < 4; i++)
        dblock[i].Log = (Edge_*)malloc(MAX_LOG * sizeof(Edge_));

    // 多线程读数据
    pthread_t thread[4];
    for (int i = 0; i < 4; i++)  pthread_create(&thread[i], NULL, read, dblock + i);
    for (int i = 0; i < 4; i++)  pthread_join(thread[i], NULL);
    max_edge = dblock[0].log_len+dblock[1].log_len+dblock[2].log_len+dblock[3].log_len;
    printf("log len: %d\n", max_edge);
    bool flag = 1;
    for (int i = 0; i < 4; i++)
    {
        Edge_* Log = dblock[i].Log;
        for (int j = 0; j < dblock[i].log_len; j++)
        {
            Edge_ e = Log[j];
            if (e.u >= Max_ID || e.v >= Max_ID){
                flag = 0;
                break;
            }
        }
        if(!flag) break;
    }

    // 判断是否可以用数组做映射
    if(flag){
        // 数组映射
        Edge_ e;
        for (int i = 0; i < 4; i++)
        {
            Edge_* Log = dblock[i].Log;
            for (int j = 0; j < dblock[i].log_len; j++)
            {
                e = Log[j];
                id2num[e.u] = 1, id2num[e.v] = 1;
            }
        }
        for(int i = 0; i < Max_ID; i++)
            if (id2num[i]) l1[max_id] = i, id2num[i] = max_id++;

        // 对结点hash
        for (int i = 0; i < 4; i++)  pthread_create(&thread[i], NULL, arrayHash, dblock + i);
        for (int i = 0; i < 4; i++)  pthread_join(thread[i], NULL);
    }
    else{
        // map映射
        Edge_ e;
        for (int i = 0; i < 4; i++)
        {
            Edge_* Log = dblock[i].Log;
            int len = dblock[i].log_len;
            for (int j = 0; j < len; j++)
            {
                Edge_ e = Log[j];
                int u = e.u, v = e.v;
                if(id2numMap.find(u) == id2numMap.end()) id2numMap[u] = 1, l1[max_id++] = u;
                if(id2numMap.find(v) == id2numMap.end()) id2numMap[v] = 1, l1[max_id++] = v;
            }
        }
        sort(l1, l1 + max_id);
        for(int i = 0; i < max_id; i++) id2numMap[l1[i]] = i;
        // 对结点hash
        for (int i = 0; i < 4; i++)  pthread_create(&thread[i], NULL, mapHash, dblock + i);
        for (int i = 0; i < 4; i++)  pthread_join(thread[i], NULL);
    }
    printf("max id: %d\n", max_id);

    pthread_create(&thread[0], NULL, get_strmap, NULL);
    pthread_create(&thread[1], NULL, create_graph, NULL);
    pthread_join(thread[0], NULL);
    pthread_join(thread[1], NULL);
}

inline bool checkw(LL a, LL b)
{
    return 5 * b < a || b > 3 * a;
}

void insert_sort(Node a[],  int length)
{
	for (int i = 1; i < length; i++)
	{
		for (int j = i - 1; j >= 0 && a[j + 1].n2 < a[j].n2; j--)
		{
            Node tmp = a[j];
            a[j] = a[j + 1];
            a[j + 1] = tmp;
		}
	}
}

void *find(void* arg)
{
    /*
        找环， 先反向找三层，再正向找四层
        为了提高取址速度，对菊花点和非菊花点动态切换两种存储方式
        为了提高程序效率，将dfs转为迭代
        由于path数组有序程度较高，使用插入排序来提高排序效率
        path数组映射版慢于不映射版，所以A榜时没有使用映射
        为了提高效率，尽可能压缩数据类型的字节数
    */
    Timer timer;
    CBlock* blk = (CBlock*)arg;
    bool *mark = (bool*)calloc(max_id, sizeof(bool));
    Node *path = (Node*)malloc(max_id * path_size * sizeof(Node));
    uint8_t *n_path = (uint8_t*)calloc(max_id, sizeof(uint8_t));
    int *history = (int*)malloc(max_id * sizeof(int));
    int n_hist = 0;
    int n3 = 0, n4 = 0, n5 = 0, n6 = 0, n7 = 0;
    int *cc3 = blk->ccs[0], *cc4 = blk->ccs[1], *cc5 = blk->ccs[2], *cc6 = blk->ccs[3], *cc7 = blk->ccs[4];
    int *nc3 = blk->ncs[0], *nc4 = blk->ncs[1], *nc5 = blk->ncs[2], *nc6 = blk->ncs[3], *nc7 = blk->ncs[4];
    int last_id1 = 0;
    int id1, id2, id3, id4, id5, id6, id7;
    Edge *j, *j2, *j3, *j4, *j5, *j6, *j7, edge_tmp;
    LL w2, w3, w4, w5, w6, w7, w8;
    Node *tmp;
    US tar2, tar3, tar4, tar5;
    while (true)
    {
        pthread_mutex_lock(&mutex);
        id1 = task;
        task += block_size;
        pthread_mutex_unlock(&mutex);
        if (id1 >= max_id) break;
        fill(mark + last_id1, mark + id1, 1);
        last_id1 = min(id1 + block_size - 1, max_id);
        for (int t = 0; t < block_size && id1 < max_id; t++, id1++)
        {
            blk->ids[blk->n_id++] = id1;
            mark[id1] = 1;

            /// 反向找三层
            for (int i = 0; i < n_hist; i++) n_path[history[i]] = 0;
            n_hist = 0;
            if (isDenseGb[id1]){
                j2 = denseGb + gbidx[id1];
                tar2 = denseIdeg[id1];
            }
            else{
                j2 = Gb + id1 * T;
                tar2 = ideg[id1];
            }
            for (int i2 = 0; i2 < tar2; i2++)
            {
                edge_tmp = j2[i2];
                id2 = edge_tmp.v;
                if (mark[id2]) continue;
                w2 = edge_tmp.w;
                mark[id2] = 1;
                if (isDenseGb[id2]){
                    j3 = denseGb + gbidx[id2];
                    tar3 = denseIdeg[id2];
                }
                else{
                    j3 = Gb + id2 * T;
                    tar3 = ideg[id2];
                }
                for (int i3 = 0; i3 < tar3; i3++)
                {
                    edge_tmp = j3[i3];
                    id3 = edge_tmp.v;
                    if (mark[id3]) continue;
                    w3 = edge_tmp.w;
                    if (checkw(w3, w2)) continue;
                    if (isDenseGb[id3]){
                        j4 = denseGb + gbidx[id3];
                        tar4 = denseIdeg[id3];
                    }
                    else{
                        j4 = Gb + id3 * T;
                        tar4 = ideg[id3];
                    }
                    for (int i4 = 0; i4 < tar4; i4++)
                    {
                        edge_tmp = j4[i4];
                        id4 = edge_tmp.v;
                        if (mark[id4]) continue;
                        w4 = edge_tmp.w;
                        if (checkw(w4, w3)) continue;
                        if(n_path[id4] == 0) history[n_hist++] = id4;
                        path[id4 * path_size + n_path[id4]++] = {id2, id3, (int)w2, (int)w4};
                    }
                }
                mark[id2] = 0;
            }
            // 排序
            for(int i = 0; i < n_hist; i++){
                int id = history[i];
                insert_sort(path + id * path_size, n_path[id]);
            }

            if (isDenseG[id1]){
                j2 = denseG + gidx[id1];
                tar2 = denseOdeg[id1];
            }
            else{
                j2 = G + id1 * T;
                tar2 = odeg[id1];
            }
            // c7
            for (int i2 = 0; i2 < tar2; i2++)
            {
                id2 = j2[i2].v;
                if (mark[id2] == 1) continue;
                mark[id2] = 1;
                w2 = j2[i2].w;
                // 四元环
                tmp = path + id2 * path_size;
                for(int i = 0; i < n_path[id2]; i++){
                    w3 = tmp[i].whead, w5 = tmp[i].wtail;
                    if(!checkw(w2, w3) && !checkw(w5, w2)){
                        id3 = tmp[i].n2, id4 = tmp[i].n1;
                        int *k = cc4 + n4 * 4;
                        *(k) = id1, *(k+1) = id2, *(k+2) = id3, *(k+3) = id4;
                        n4++;  nc4[id1]++;
                    }
                }

                if (isDenseG[id2]){
                    j3 = denseG + gidx[id2];
                    tar3 = denseOdeg[id2];
                }
                else{
                    j3 = G + id2 * T;
                    tar3 = odeg[id2];
                }
                for (int i3 = 0; i3 < tar3; i3++)
                {
                    id3 = j3[i3].v, w3 = j3[i3].w;
                    if (mark[id3] == 1 || checkw(w2, w3)) continue;
                    // 五元环
                    tmp = path + id3 * path_size;
                    for(int i = 0; i < n_path[id3]; i++){
                        id4 = tmp[i].n2, id5 = tmp[i].n1;
                        if(id4 == id2 || id5 == id2 ) continue;
                        w4 = tmp[i].whead, w6 = tmp[i].wtail;
                        if(!checkw(w3, w4) && !checkw(w6, w2)){
                            int *k = cc5 + n5 * 5;
                            *(k) = id1, *(k+1) = id2, *(k+2) = id3, *(k+3) = id4, *(k+4) = id5;
                            n5++;  nc5[id1]++;
                        }
                    }
                    mark[id3] = 1;

                    if (isDenseG[id3]){
                        j4 = denseG + gidx[id3];
                        tar4 = denseOdeg[id3];
                    }
                    else{
                        j4 = G + id3 * T;
                        tar4 = odeg[id3];
                    }
                    for (int i4 = 0; i4 < tar4; i4++)
                    {
                        id4 = j4[i4].v, w4 = j4[i4].w;
                        if (checkw(w3, w4)) continue;
                        if(id4 == id1 && !checkw(w4, w2)){
                            int *k = cc3 + n3 * 3;
                            *(k) = id1;  *(k+1) = id2;  *(k+2) = id3;
                            n3++;  nc3[id1]++;
                        }
                        if (mark[id4] == 1) continue;
                        // 六元环
                        tmp = path + id4 * path_size;
                        for(int i = 0; i < n_path[id4]; i++){
                            id5 = tmp[i].n2, id6 = tmp[i].n1;
                            if(mark[id5] || mark[id6]) continue;
                            w5 = tmp[i].whead, w7 = tmp[i].wtail;
                            if(!checkw(w4, w5) && !checkw(w7, w2)){
                                int *k = cc6 + n6 * 6;
                                *(k) = id1, *(k+1) = id2, *(k+2) = id3, *(k+3) = id4, *(k+4) = id5, *(k+5) = id6;
                                n6++;  nc6[id1]++;
                            }
                        }
                        mark[id4] = 1;

                        if (isDenseG[id4]){
                            j5 = denseG + gidx[id4];
                            tar5 = denseOdeg[id4];
                        }
                        else{
                            j5 = G + id4 * T;
                            tar5 = odeg[id4];
                        }
                        for (int i5 = 0; i5 < tar5; i5++)
                        {
                            id5 = j5[i5].v, w5 = j5[i5].w;
                            if (n_path[id5] == 0 || mark[id5] || checkw(w4, w5)) continue;
                            // 七元环
                            tmp = path + id5 * path_size;
                            for(int i = 0; i < n_path[id5]; i++){
                                id6 = tmp[i].n2, id7 = tmp[i].n1;
                                if (mark[id6] || mark[id7]) continue;
                                w6 = tmp[i].whead, w8 = tmp[i].wtail;
                                if(!checkw(w5, w6) && !checkw(w8, w2)){
                                    int *k = cc7 + n7 * 7;
                                    *(k) = id1, *(k+1) = id2, *(k+2) = id3, *(k+3) = id4, *(k+4) = id5, *(k+5) = id6, *(k+6) = id7;
                                    n7++;  nc7[id1]++;
                                }
                            }
                        }
                        mark[id4] = 0;
                    }
                    mark[id3] = 0;
                }
                mark[id2] = 0;
            }
        }
    }
    blk->ns[0] = n3; blk->ns[1] = n4; blk->ns[2] = n5; blk->ns[3] = n6; blk->ns[4] = n7;
    // 确定ID在cc中的位置
    for (int i = 0; i < 5; i++)
    {
        int l = i + 3, *ic = blk->ics[i], *nc = blk->ncs[i];
        ic[0] = 0;
        for (int i = 1; i < max_id; i++) ic[i] = ic[i-1] + nc[i-1] * l;
    }
    timer.stop("find(*)");
}

void find()
{
    /*
        使用抢占式负载均衡方式多线程找环
        找完环之后，确定每个环的指针所在
    */
    // 初始化
    for (int i = 0; i < NT_F; i++)
    {
        CBlock* blk = cblock + i;
        for (int j = 0; j < 5; j++)
        {
            blk->ns[j] = 0;
            blk->ccs[j] = (int*)malloc(20000000 * (j+3) * sizeof(int));
            blk->ncs[j] = (int*)calloc(max_id, sizeof(int));
            blk->ics[j] = (int*)malloc(max_id * sizeof(int));
        }
        blk->n_id = 0;
        blk->ids = (int*)malloc(max_id * sizeof(int));
    }
    // 找环
    task = 0;
    pthread_t thread[NT_F];
    for (int i = 0; i < NT_F; i++)  pthread_create(&thread[i], NULL, find, cblock + i);
    for (int i = 0; i < NT_F; i++)  pthread_join(thread[i], NULL);
    // 确定ID在哪个block中
    uint8_t *in_blk = (uint8_t*)malloc(max_id * sizeof(uint8_t));
    memset(in_blk, 100, max_id);
    for (int i = 0; i < NT_F; i++)
    {
        CBlock* blk = cblock + i;
        for (int j = 0; j < 5; j++)
        {
            total_r += blk->ns[j];
            each_r[j] += blk->ns[j];
        }
        for (int j = 0; j < blk->n_id; j++)
        {
            in_blk[blk->ids[j]] = i;
        }
    }
    // 确定所有环路的指针
    pos = (int**)malloc(total_r * sizeof(int*));
    plen = (int*)malloc(total_r * sizeof(int));
    int ip = 0, *p, nn;
    CBlock *blk;
    for (int i = 0; i < 5; i++)
    {
        int l = i + 3;
        for (int j = 0; j < max_id; j++)
        {
            if (in_blk[j] == 100) continue;
            blk = cblock + in_blk[j];
            p = blk->ccs[i] + blk->ics[i][j];
            nn = blk->ncs[i][j];
            for (int k = 0; k < nn; k++)
            {
                plen[ip] = l;
                pos[ip++] = p;
                p += l;
            }
        }
    }
    printf("total_r: %d\neach: ", total_r);
    for (int i = 0; i < 5; i++) printf("%d ", each_r[i]);
}

void *circle2str(void* arg)
{
    /*
        将结果转为字符串并写入
        利用标记来确认是否已经获得前面环的总字节长度
        memcpy利用对齐技巧，每次都指定长度为16，下次按照真实长度开始覆盖
    */
    Timer timer;
    CRange* range = (CRange*)arg;
    int maxcsz = (total_r / 8 + 1) * 7 * 11;
    int i = range->i, begin = range->begin, end = range->end, csize = 0;
    char* str = (char*)malloc(maxcsz * sizeof(char));
    for (int j = begin; j < end; j++)
    {
        int *p = pos[j], l = plen[j];
        for (int k = 0; k < l; k++)
        {
            int c = p[k];
            char* ch = node2str + c * 16;
            int lc = node2len[c];
            memcpy(str+csize, ch, 16);
            csize += lc;
        }
        str[csize-1] = '\n';
    }
    timer.stop("c2s*");
    int wp;
    if (i == 0) {
        wp = fsize; fsize += csize; done[0] = 1;
        pwrite(fd, str, csize, wp);
    } else {
        while (1) if (done[i-1]) {wp = fsize; fsize += csize; done[i] = 1; break;}
        pwrite(fd, str, csize, wp);
    }
}

void save()
{
    /*
        将结果转为字符串并写入
    */
    // 总数
    char* odata = (char*)malloc(10 * sizeof(char));
    sprintf(odata, "%d\n", total_r);
    fsize = strlen(odata);
    fd = open(output_path, O_RDWR | O_CREAT, 0666);
    write(fd, odata, fsize);
    // 按比例均分并转字符串+写入
    float ratio[8] = {0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
    crange[0].begin = 0;
    for (int i = 0; i < 8; i++)
    {
        crange[i].i = i;
        crange[i].end = int(ratio[i] * total_r);
        done[i] = 0;
    }
    for (int i = 1; i < 8; i++)  crange[i].begin = crange[i-1].end;
    crange[7].end = total_r;
    pthread_t thread[8];
    for (int i = 0; i < 8; i++) pthread_create(&thread[i], NULL, circle2str, crange + i);
    for (int i = 0; i < 8; i++)  pthread_join(thread[i], NULL);
    printf("fsize: %d\n", fsize);
}

int main(int argc, char *argv[])
{
    Timer timer1;
    load();
    timer1.stop("load");

    Timer timer2;
    find();
    timer2.stop("find");

    Timer timer3;
    save();
    timer3.stop("save");
    return 0;
}
