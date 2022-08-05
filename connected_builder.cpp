#include <bits/stdc++.h>
#include "Eigen/Eigenvalues"
using namespace std;
 
typedef long long ll;
typedef vector<int> vi;
typedef vector<ll> vll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
 
#define FOR(i, a, b) for (int i=a; i<(b); i++)
#define F0R(i, a) for (int i=0; i<(a); i++)
#define FORd(i,a,b) for (int i = b; i >= a; i--)
#define F0Rd(i,a) for (int i = a; i >= 0; i--)
#define FORit(it,a) for (auto it = a.begin(); it != a.end(); it++)
#define trav(a,x) for (auto& a: x)
 
#define sz(x) (int)(x).size()
#define mp make_pair
#define pb push_back
#define all(x) x.begin(), x.end()
#define alla(arr, sz) arr, arr + sz
 
const int dx[4] = { 1,0,-1,0 }, dy[4] = { 0,1,0,-1 }; // for every grid problem!!
const ll linf = 4*1e18;
const ll inf = 1000000007;
 
namespace io {
	void setIn(string s) { freopen(s.c_str(),"r",stdin); }
	void setOut(string s) { freopen(s.c_str(),"w",stdout); }
	void setIO(string s = "") {
		ios_base::sync_with_stdio(0); cin.tie(0); // fast I/O
		if (sz(s)) { setIn(s+".in"), setOut(s+".out"); } // for USACO
	}
}
 
using namespace io;

#define STRONG 1

const int mxn = 1e5+5;

int n, root; 
int global_time;

int startTime[mxn];
int low[mxn];
int inverse_startTime[mxn];

vector<vi> adj;
int edge_cnt;

vi st;
void build_st() {
    st.pb(root);
    F0R(i,n) if (i != root) {
        int l = st.size();
        int ridx = rand() % l;
        adj[st[ridx]].pb(i);
        st.pb(i);
    }
}

void tarjan(int v) {
    startTime[v] = global_time;
    low[v] = startTime[v];
    inverse_startTime[global_time] = v;
    global_time++;

    trav(x,adj[v]) {
        if (startTime[x] == -1) {
            tarjan(x);
            low[v] = min(low[v], low[x]);
        }
        low[v] = min(low[v], startTime[x]);
    }

    if (low[v] == startTime[v] && low[v] != 0) {
        int a = rand() % startTime[v]; // smaller start time: [0, startTime[v])
        int b = rand() % (global_time - startTime[v]) + startTime[v]; // same or greater start time: [startTime[v], global_time)
        int x = inverse_startTime[a];
        int y = inverse_startTime[b];
        adj[y].pb(x);
        low[v] = a;
        edge_cnt++;
    }
}

vector<pii> edges;

void shuffle() {
    for (int i = 0; i < edges.size()-1; i++) {
        int j = rand() % (edges.size()-i) + i; // [i, edges.size())
        auto t = edges[i];
        edges[i] = edges[j];
        edges[j] = t;
    }
}

void erdos_renyi(int m) {

    F0R(i,n) F0R(j,n) if (i != j) {
        bool b = true;
        trav(x, adj[i]) {
            if (x == j) {
                b = false; continue;
            }
        }
        if (b) edges.pb(mp(i,j));
    }

    shuffle();
    
    F0R(i,m) {
        pii t = edges[i];
        adj[t.first].pb(t.second);
    }
}

void add_edges(){
    F0R(i,n) F0R(j,n) if (i != j) {
        bool b = true;
        trav(x, adj[i]) {
            if (x == j) {
                b = false; continue;
            }
        }
        if (b) {
            int p = rand() % 2; // 0 to 1
            if (p == 1) {
                adj[i].pb(j);
                edge_cnt++;
            }
        }
    }
}

void init() {
    global_time = 0;
    memset(startTime, -1, sizeof startTime);
    memset(low, 0, sizeof low);
    memset(inverse_startTime, 0, sizeof inverse_startTime);
    adj.clear();
    edge_cnt = 0;
    st.clear();
    edges.clear();
}

int main() {
    #if STRONG
    //setOut("strong.txt");
    #else
    setOut("root.txt");
    #endif

    srand (time(NULL));

    FOR(node_cnt,5,10) { F0R(iteration,50) {
        init();

        n = node_cnt;
        // n = 5;

        adj.resize(n);
        //root = rand() % n;
        root = 0;
        // printf("root: %d\n", root);
        build_st();

        edge_cnt = n-1;

        #if STRONG
        tarjan(root); // tarjan function updates edge_cnt
        #endif

        int additional_edges = rand() % (n*(n-1)-edge_cnt+1); // [0, n*(n-1)-edge_cnt]
        erdos_renyi(additional_edges); 

        // F0R(i,n) {
        //     cout << i << ": "; 
        //     trav(x,adj[i]) {
        //         cout << x << ' ';
        //     }
        //     cout << '\n';
        // }
        
        edge_cnt += additional_edges; // add the additional edges, now edge_cnt should be finalized 
        Eigen::MatrixXd B(n, edge_cnt);
        
        int flag = 0; // column index in B
        F0R(i,n) {
            trav(x,adj[i]) {
                F0R(k,n) B(k,flag) = 0; // set entire column to 0

                B(i,flag) = 1;
                B(x,flag) = -1;
                flag++; // move to next column
            }
        }

        Eigen::MatrixXd B_T = B.transpose();
        Eigen::MatrixXd L = B*B_T;

        Eigen::EigenSolver<Eigen::MatrixXd> es;
        es.compute(L, /* computeEigenvectors = */ false);
        auto e = es.eigenvalues(); // returns a column vector of type complex<double>

        vector<double> eigenval;
        trav(x, e) eigenval.pb(x.real());
        sort(all(eigenval));
        cout << eigenval[1] << ' '; 

        /*
          add 1 more edge

        */

        B.conservativeResize(B.rows(), edge_cnt+1);

        FOR(i, additional_edges, edges.size()) { // edges.size() = n*(n-1)-edge_cnt
            F0R(j,n) B(j,flag) = 0;
            pii next = edges[i];
            B(next.first, flag) = 1;
            B(next.second, flag) = -1;

            B_T = B.transpose();
            L = B*B_T;
            es.compute(L, false);
            e = es.eigenvalues();
            vector<double> eigenval;
            trav(x, e) eigenval.pb(x.real());
            sort(all(eigenval));
            cout << eigenval[1] << ' '; 
        }

        cout << endl;
        
    } cout << endl; }

    return 0;
}