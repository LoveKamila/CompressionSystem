#include <fstream>
#include <iostream>
#include<cstdio>
#include <cstdlib>
#include <string>
#include<vector>
#include<queue>
#include<map>
#include <bitset>
using namespace std;

//Huffman���Ľ��
struct node {
    long long freq;  //�ַ�����Ƶ��
    unsigned char chr;   //�ַ�
    node* left;  //����
    node* right;  //�Һ���
    node* father;  //���ڵ�
    node() {
        chr = 0;
        left = NULL;
        right = NULL;
        father = NULL;
        freq = -1;
    }
};

//���ȶ���ά����㣬�����ַ�Ƶ����С��������
struct cmp {
    bool operator()(const node& n1, const node& n2) {
        return n1.freq > n2.freq;
    }
};
priority_queue <node, vector<node>, cmp> que;

long long cnt[256] = { 0 };  //ͳ��ÿ���ַ�Ƶ��
string content;  //����ѹ�����ļ�һ���Զ����ַ�����

struct CODE {
    int len = 0;  //���ַ��ı��볤��
    bitset<32> cd = 0;  //���ַ��ı���
};
struct CODE code[256];  //���ַ���ascii��Ϊ�±꣬�������


//����Huffman��
struct node* build_tree()
{
    //�����ֹ����ַ��������ȶ���
    for (int i = 0; i < 256; i++)
        if (cnt[i] > 0) {
            struct node* newnode = new(nothrow)(struct node);
            if (newnode == NULL) {
                cout << "wrong" << endl;
                exit(-1);
            }
            newnode->freq = cnt[i];
            newnode->chr = (unsigned char)i;
            que.push(*newnode);
        }
    //������������
    while (!que.empty()) {
        struct node* top1 = new(nothrow)(struct node);
        struct node* top2 = new(nothrow)(struct node);
        struct node* fa = new(nothrow)(struct node);  //���ڵ�
        if (top1 == NULL || top2 == NULL || fa == NULL) {
            cout << "wrong" << endl;
            exit(-1);
        }
        *top1 = que.top();
        que.pop();
        if (que.empty())
            return top1;
        *top2 = que.top();
        que.pop();
        top1->father = fa;
        top2->father = fa;
        fa->left = top1;
        fa->right = top2;  
        fa->freq = top1->freq + top2->freq;
        que.push(*fa);
    }
    return NULL;
}

//�õ����ַ��ı���
void get_code
(struct node* root, int len, bitset<32> x)
{
    //����Ҷ�ӽ��
    if (root->left == NULL && root->right == NULL) {
        int i = (int)root->chr;
        code[i].len = len;
        code[i].cd = x;
        return;
    }
    //����ݹ�
    x[len] = 0;
    get_code(root->left, len + 1, x);
    //���ҵݹ�
    x[len] = 1;
    get_code(root->right, len + 1, x);
    return;
}

//����Huffman���ռ�
void cleartree(node* root)
{
    if (!root)
        return;
    cleartree(root->left);
    cleartree(root->right);
    delete(root);
}


int main(int argc, char* argv[]) {
    //���������������ʽ�ж�
    cout << "Zipper 0.001! Author: root" << endl;
    if (argc != 4) {
        cerr << "Please make sure the number of parameters is correct." << endl;
        return -1;
    }
    if (strcmp(argv[3], "zip") && strcmp(argv[3], "unzip")) {
        cerr << "Unknown parameter!\nCommand list:\nzip/unzip" << endl;
        return -1;
    }

    //�ļ���
    ifstream fin(argv[1], ios::binary); // �Զ����Ʒ�ʽ�������ļ�
    if (!fin) {  // ���������Ϣ���˳�
        cerr << "Can not open the input file!" << endl;
        return -1;
    }
    ofstream fout(argv[2], ios::binary); // �Զ����Ʒ�ʽ������ļ�
    if (!fout) {  // ���������Ϣ���˳�
        cerr << "Can not open the output file!" << endl;
        return -1;
    }

    //����ѹ�����ѹ��
    /*****************����ѹ������******************/
    if (strcmp(argv[3], "zip") == 0) {
        istreambuf_iterator<char> beg(fin), end; 
        // ���������ļ�ָ�룬ָ��ʼ�ͽ������� char(һ�ֽ�) Ϊ����
        string content(beg, end); 
        // ���ļ�ȫ������ string �ַ���

        //ͳ�Ƹ��ַ�Ƶ��
        int pos = 0;
        while (pos < (int)content.length()) {
            cnt[(int)content[pos]]++;
            pos++;
        }

        //���������
        bitset<32> x = 0;  //����ÿ���ַ���huffman����
        struct node* root = build_tree();   //����Huffman��
        get_code(root, 0, x);   //�õ����ַ�����
        
        //��8��������λΪ��λд���ļ�
        bitset<8> n = 0;  
        int bitnum = 0;  //��ǰ������λ����bitset��λ��
        pos = 0;  //ԭ�ļ��ַ����±�
        while (pos < (int)content.length())
        {
            int asc = (int)content[pos];
            int i = code[asc].len;  //��ǰ�ַ���01�����볤��
            for (int j = 0; j < i; j++) {
                n[bitnum] = code[asc].cd[j];
                bitnum++;
                if (bitnum >= 8) {  //д��8λ�����
                    fout.put((char)n.to_ulong());
                    n = 0;
                    bitnum = 0;
                }
            }
            pos++;
        }
        if (bitnum != 0)   //������ʣ�಻��8λ��
            fout.put((char)n.to_ulong());

        //������ַ�Ƶ�����ڽ�ѹ��ʱ��ԭ
        ofstream outtree("tree.log", ios::out);
        if (!outtree) {
            cerr << "Can not open the output file!" << endl;
            return -1;
        }
        for (int i = 0; i < 256; i++)
            outtree << cnt[i] << " ";
        outtree.close();

        //�ͷ�Huffman���ռ�
        cleartree(root);
    }

    /*****************��ѹ������******************/
    else if (strcmp(argv[3], "unzip") == 0) {
        //������ַ���Ƶ��
        ifstream intree("tree.log", ios::in);
        if (!intree) {
            cerr << "Can not open the output file!" << endl;
            return -1;
        }
        for (int i = 0; i < 256; i++)
            intree >> cnt[i];
        intree.close();

        //���½���Huffman��
        struct node* root = build_tree();
        if (root == NULL) {
            cerr << "wrong" << endl;
            exit(-1);
        }

        //����ѹ������ļ�����ԭ����
        istreambuf_iterator<char> beg(fin), end; // ���������ļ�ָ�룬ָ��ʼ�ͽ������� char(һ�ֽ�) Ϊ����
        string s(beg, end); // ���ļ�ȫ������ string �ַ���
        int pos = 0;
        struct node* tmproot = root;  //��ʾ��ǰ�ڵ�
        while (pos < (int)s.length()) {
            bitset<8> x((unsigned long)s[pos]); //����ǰ�ַ�תΪ01��
            int i = 0;
            while (i < 8) {
                //��0�������ӣ�1�����Һ���
                if (x[i] == 0)
                    tmproot = tmproot->left;
                else
                    tmproot = tmproot->right;
                i++;
                //����Ҷ�ӽ��������·������ʾ���ַ�
                if (tmproot->left == NULL && tmproot->right == NULL) {
                    fout.put(tmproot->chr);
                    tmproot = root;
                }
            }
            pos++;
        }

        //�ͷ�Huffman���ռ�
        cleartree(root);
    }

    // ȫ���������ļ���ر��ļ����
    fin.close();
    fout.close();

    cout << "Complete!" << endl;

    return 0;
}