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

using namespace std;

#define MAX_DEG 32
#define MAX_ID 50001                        // 根据大佬提示，线上有效ID不超过50000
// #define MAX_ID 200000
#define NT_F 8  // 找环线程数

const char* input_path = "/data/test_data.txt";
const char* output_path = "/projects/student/result.txt";

int max_id;                                 // 最大ID
char** node2str;                            // ID映射字符串
char* node2len;                             // ID映射字符串长度
int *G, *Gb;                                // 正、反向邻接表
uint8_t *ideg, *odeg;                       // 出入度

// 找环数据块：记录单个线程找环所用变量
struct CBlock
{
    int begin, end;                         // 找环的起点、终点ID
    int n3, n4, n5, n6, n7;                 // 每种环的数量
    int *cc3, *cc4, *cc5, *cc6, *cc7;       // 保存每种环的结果
    int *ncs;                               // {n3, n4, n5, n6, n7}
    int **ccs;                              // {*cc3, *cc4, *cc5, *cc6, *cc7}
};
CBlock cblock[NT_F];

// 转换字符串范围：记录单个线程转换字符串的数据范围
struct CRange  
{
    int i;                                  // 序号
    int g;                                  // 从第g个的数组开始
    int c;                                  // 从该数组第c个环开始
    int l;                                  // 转换的环数
};
CRange crange[8];

volatile uint8_t done[8];                   // 是否完成转换字符串的标志
volatile int fsize;                         // 结果文件大小
int fd;                         

// 初始化变量
void init()
{
    node2str = (char**)malloc(MAX_ID * sizeof(char*));
    node2len = (char*)malloc(MAX_ID * sizeof(char));
    G = (int*)malloc(MAX_ID * MAX_DEG * sizeof(int));
    Gb = (int*)malloc(MAX_ID * MAX_DEG * sizeof(int));
    ideg = (uint8_t*)calloc(MAX_ID, sizeof(uint8_t));
    odeg = (uint8_t*)calloc(MAX_ID, sizeof(uint8_t));
    // b
    for (int i = 0; i < NT_F; i++)
    {   
        CBlock* blk = cblock + i;
        blk->n3 = blk->n4 = blk->n5 = blk->n6 = blk->n7 = 0;
        blk->cc3 = (int*)malloc(750000 * sizeof(int)); // 500000 * 3 / 2
        blk->cc4 = (int*)malloc(1000000 * sizeof(int)); // 500000 * 4 / 2
        blk->cc5 = (int*)malloc(2500000 * sizeof(int)); // 1000000 * 5 / 2
        blk->cc6 = (int*)malloc(6000000 * sizeof(int)); // 2000000 * 6 / 2
        blk->cc7 = (int*)malloc(10500000 * sizeof(int)); // 3000000 * 7 / 2
    }
}

// 整数转字符串
void itoa(int d, char** c, char* l)
{
    char* s = (char*)malloc(10*sizeof(char));
    char i = 9;
    s[i--] = ',';
    while (d > 0)
    {
        s[i--] = d % 10 + 48;
        d /= 10;
    }
    *c = s + i + 1;
    *l = 9 - i;
}

// 获取文件大小
int get_file_size(const char *path)  
{       
    struct stat statbuff; 
    stat(path, &statbuff); 
    return statbuff.st_size;  
}

// 映射字符串
void *get_strmap(void* arg)
{
    node2str[0] = (char*)calloc(1, sizeof(char)); 
    node2str[0][0] = '0'; node2str[0][1] = ','; node2len[0] = 2;
    for (int i = 1; i < MAX_ID; i++)  itoa(i, node2str + i, node2len + i);
}

// 加载数据、建图
void* load(void* arg)
{
    // 读取
    int ifsize = get_file_size(input_path);
    int fd = open(input_path, O_RDONLY, S_IRUSR);
    char *ifmap = (char *)mmap(NULL, ifsize, PROT_READ, MAP_SHARED, fd, 0);
    char c;
    int i = 0, node1, node2;   
    while (true)
    {
        node1 = 0; node2 = 0;
        while ((c = ifmap[i++]) != ',')  node1 = (node1 * 10 + c - 48);
        while ((c = ifmap[i++]) != ',')  node2 = (node2 * 10 + c - 48);
        while ((c = ifmap[i++]) != '\n');
        if (node1 < MAX_ID && node2 < MAX_ID)  
        {
            G[node1 * MAX_DEG + odeg[node1]++] = node2;
            Gb[node2 * MAX_DEG + ideg[node2]++] = node1;
        }
        if (i >= ifsize) break;
    }
    // 确定实际最大id
    for (i = MAX_ID-1; i >=0; i--) {if(odeg[i] != 0) break;}
    max_id = i + 1;
    printf("max_id: %d\n", max_id);
    // 排序
    int j, k;
    int *p1, *p2;
    for (i = 0; i < max_id; i++)
    {
        p1 = G + i * MAX_DEG; p2 = Gb + i * MAX_DEG;
        sort(p1, p1 + odeg[i]);
        sort(p2, p2 + ideg[i]);
    }
    // 确定找环区间，按ID均匀分配
    float ratio[NT_F] = {0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1}; 
    for (i = 0; i < NT_F; i++)
    {
        cblock[i].begin = 0;
        cblock[i].end = int(ratio[i] * max_id);
    }
    for (int i = 1; i < NT_F; i++)  cblock[i].begin = cblock[i-1].end;
}

// 多线程：一个线程加载数据、建图，一个线程映射字符串
void load()
{
    pthread_t thread[2];
    pthread_create(&thread[1], NULL, load, NULL);
    pthread_create(&thread[0], NULL, get_strmap, NULL);
    for (int i = 0; i < 2; i++)  pthread_join(thread[i], NULL);
}

// BFS反向遍历三层，记录节点到起点的距离（深度）
void get_depth(int node, int *queue, uint8_t *depth, int *history, int &n_hist)
{
    for (int i = 0; i < n_hist; i++) depth[history[i]] = 10; 
    depth[node] = 0; queue[0] = node;
    int d = 1, p1 = 0, p2 = 0, p3 = 1;
    n_hist = 0;
    int *j;
    while (p1 < p3)
    {
        p2 = p3;
        while (p1 < p2)
        {
            int curr = queue[p1++];
            j = Gb + curr * MAX_DEG;
            for (int i = 0; i < ideg[curr]; i++)
            {
                int next = *(j + i);
                if (depth[next] == 10) 
                {
                    depth[next] = d;
                    queue[p3++] = next;
                    history[n_hist++] = next;
                }   
            }
        }
        d++;
        if (d > 3) break;
    }
    depth[node] = 9;
}

// 找环：先反向遍历三层记录深度，再前向遍历6层，利用深度进行减枝
void *find(void* arg)
{
    CBlock* blk = (CBlock*)arg;
    int begin = blk->begin, end = blk->end;
    uint8_t *mark = (uint8_t*)calloc(max_id, sizeof(uint8_t));
    int *queue = (int*)malloc(max_id * sizeof(int));
    int *history = (int*)malloc(max_id * sizeof(int));
    uint8_t *depth = (uint8_t*)malloc(max_id * sizeof(uint8_t));
    memset(depth, 10, max_id);
    memset(depth, 9, begin);
    memset(mark, 1, begin);
    int n_hist = 0;
    for (int node1 = begin; node1 < end; node1++)
    {
        get_depth(node1, queue, depth, history, n_hist);
        mark[node1] = 1;
        int* j2 = G + node1 * MAX_DEG;  
        // c7
        for (int i2 = 0; i2 < odeg[node1]; i2++)
        {
            int node2 = *(j2 + i2);
            if (mark[node2] == 1) continue;
            mark[node2] = 1;
            int* j3 = G + node2 * MAX_DEG;
            for (int i3 = 0; i3 < odeg[node2]; i3++)
            {
                int node3 = *(j3 + i3);
                if (mark[node3] == 1) continue;
                if (depth[node3] == 1)
                {
                    int *k = blk->cc3 + blk->n3 * 3;
                    *(k) = node1;  *(k+1) = node2;  *(k+2) = node3;
                    blk->n3++;
                }
                mark[node3] = 1;
                int* j4 = G + node3 * MAX_DEG;
                for (int i4 = 0; i4 < odeg[node3]; i4++)
                {
                    int node4 = *(j4 + i4);
                    if (mark[node4] == 1) continue;
                    if (depth[node4] == 1)
                    {
                        int *k = blk->cc4 + blk->n4 * 4;
                        *(k) = node1;  *(k+1) = node2;  *(k+2) = node3; *(k+3) = node4;
                        blk->n4++;
                    }
                    mark[node4] = 1;
                    int* j5 = G + node4 * MAX_DEG;
                    for (int i5 = 0; i5 < odeg[node4]; i5++)
                    {
                        int node5 = *(j5 + i5);
                        if (depth[node5] > 3 || mark[node5] == 1) continue;
                        if (depth[node5] == 1)
                        {
                            int *k = blk->cc5 + blk->n5 * 5;
                            *(k) = node1;  *(k+1) = node2;  *(k+2) = node3; *(k+3) = node4; *(k+4) = node5;
                            blk->n5++;
                        }
                        mark[node5] = 1;
                        int* j6 = G + node5 * MAX_DEG;
                        for (int i6 = 0; i6 < odeg[node5]; i6++)
                        {
                            int node6 = *(j6 + i6);
                            if (depth[node6] > 2 || mark[node6] == 1) continue;
                            if (depth[node6] == 1)
                            {
                                int *k = blk->cc6 + blk->n6 * 6;
                                *(k) = node1;  *(k+1) = node2;  *(k+2) = node3; *(k+3) = node4; *(k+4) = node5; *(k+5) = node6;
                                blk->n6++;
                            }
                            mark[node6] = 1;
                            int* j7 = G + node6 * MAX_DEG;
                            for (int i7 = 0; i7 < odeg[node6]; i7++)
                            {
                                int node7 = *(j7 + i7);
                                if (depth[node7] > 1 || mark[node7] == 1) continue;
                                int *k = blk->cc7 + blk->n7 * 7;
                                *(k) = node1;  *(k+1) = node2;  *(k+2) = node3; *(k+3) = node4; *(k+4) = node5; *(k+5) = node6; *(k+6) = node7;
                                blk->n7++;
                            }
                            mark[node6] = 0;
                        }
                        mark[node5] = 0;
                    }
                    mark[node4] = 0;
                }
                mark[node3] = 0;
            }
            mark[node2] = 0;
        } 
    }
}

// 开多线程找环
void find()
{
    pthread_t thread[NT_F];
    for (int i = 0; i < NT_F; i++)  pthread_create(&thread[i], NULL, find, cblock + i);
    for (int i = 0; i < NT_F; i++)  pthread_join(thread[i], NULL);
    for (int i = 0; i < NT_F; i++)
    {
        CBlock* blk = cblock + i;
        blk->ncs = new int[5]{blk->n3, blk->n4, blk->n5, blk->n6, blk->n7};
        blk->ccs = new int*[5]{blk->cc3, blk->cc4, blk->cc5, blk->cc6, blk->cc7};
    }
}

// 确定每个线程转换字符串的范围，按照数字个数均匀分配
void get_range()
{
    for (int i = 0; i < 8; i++)
    {
        crange[i].i = i; crange[i].g = 0; crange[i].c = 0; crange[i].l = 0;
    }
    int ns[5 * NT_F];
    int nc = 0;
    int total_num = 0;
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < NT_F; j++)
        {
            ns[i * NT_F + j] = cblock[j].ncs[i];
            nc += cblock[j].ncs[i];
            total_num += cblock[j].ncs[i] * (i+3);
        }
    }
    float ratio[9] = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1};
    int dnum[9];
    for (int i = 0; i < 9; i++)  dnum[i] = int(total_num * ratio[i]);
    int acc_len[5 * NT_F + 1] = {0}, acc_num[5 * NT_F + 1] = {0};
    for (int i = 0; i < 5 * NT_F; i++)
    {
        acc_len[i+1] = acc_len[i] + (i / NT_F + 3)*ns[i];
        acc_num[i+1] = acc_num[i] + ns[i];        
    } 
    for (int i = 1; i < 8; i++)
    {
        int k;
        for (k = 0; k < 5 * NT_F; k++)
        {
            if (dnum[i] >= acc_len[k] && dnum[i] < acc_len[k+1]) break;
        }
        crange[i].g = k;
        crange[i].c = (dnum[i] - acc_len[k]) / (k / NT_F + 3);
        crange[i].l = crange[i].c + acc_num[k];
    }
    for (int i = 0; i < 7; i++) crange[i].l = crange[i+1].l - crange[i].l; 
    crange[7].l = nc - crange[7].l;
}

// 转换字符串并保存结果
void *circle2str(void* arg)
{
    CRange* range = (CRange*)arg;

    int ns[5 * NT_F];
    int *ccs[5 * NT_F];
    int maxcsz = 0;
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < NT_F; j++)
        {
            ns[i * NT_F + j] = cblock[j].ncs[i];
            ccs[i * NT_F + j] = cblock[j].ccs[i];
            maxcsz += cblock[j].ncs[i] * (i+3);
        }
    }
    maxcsz *= 3;
    int i = range->i, rg_g = range->g, rg_c = range->c, rg_l = range->l; 
    char* str = (char*)malloc(maxcsz * sizeof(char));
    int cnt = 0, csize = 0;
    for (int l = rg_g; l < 5 * NT_F; l++)
    {
        int len = l / NT_F + 3;
        int j = 0, mj = ns[l];
        int *r = ccs[l];
        if (l == rg_g) {j = rg_c; r += (rg_c * len);}
        while (j < mj)
        {
            for (int k = 0; k < len; k++)
            {
                int c = *r++;
                char* ch = node2str[c];
                int lc = node2len[c];
                memcpy(str+csize, ch, 8);
                csize += lc;
            }
            str[csize-1] = '\n';
            j++; cnt++;
            if (cnt >= rg_l) {l=100; break;}
        } 
    }
    int wp;
    if (i == 0) {
        wp = fsize; fsize += csize; done[0] = 1;
        pwrite(fd, str, csize, wp);
    } else {
        while (1) if (done[i-1]) {wp = fsize; fsize += csize; done[i] = 1; break;}
        pwrite(fd, str, csize, wp);
    }
}

// 开多线程转换保存
void save()
{
    // 总数
    int nc = 0;
    for (int i = 0; i < 5; i++) for (int j = 0; j < NT_F; j++) nc += cblock[j].ncs[i];
    char* odata = (char*)malloc(10 * sizeof(char));
    printf("nc: %d\n", nc);
    sprintf(odata, "%d\n", nc);
    fsize = strlen(odata);
    fd = open(output_path, O_RDWR | O_CREAT, 0666);
    write(fd, odata, fsize);
    // 转换保存
    get_range();
    done[0] = done[1] = done[2] = done[3] = done[4] = done[5] = done[6] = done[7] = 0;
    pthread_t thread[8];
    for (int i = 0; i < 8; i++) pthread_create(&thread[i], NULL, circle2str, crange + i);
    for (int i = 0; i < 8; i++)  pthread_join(thread[i], NULL);
    printf("fsize: %d\n", fsize);
}


int main(int argc, char *argv[])
{
    init();
    load();
    find();
    save();
    return 0;
}