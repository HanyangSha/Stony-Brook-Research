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

const double eps = 1e-10;
bool global_check = true; 

const int mxn = 10;
int global_time;
int startTime[mxn], low[mxn], inverse_startTime[mxn];
vector<vi> adj;
bool done;
int scc_cnt; 
vi stk; 

void init() {
    global_time = 0;
    memset(startTime, -1, sizeof startTime);
    memset(low, 0, sizeof low);
    memset(inverse_startTime, 0, sizeof inverse_startTime);
    adj.clear();
    done = false;
    scc_cnt = 0;
    stk.clear();
}

// both edge functions are 1 indexed
int add_edge(int a, int b, int n) {
    return (a-1)*n + b;
}

pii get_edge(int loc, int n) {
    int a = (int)(ceil((double)loc/n));
    int b = loc % n;
    if (loc % n == 0) b = n;

    return mp(a,b);
}

void check_rowsum(int n, Eigen::MatrixXd & B) {
    F0R(i,n) {
        int sum = 0;
        F0R(j,n) {
            sum += B(i,j);
        }
        if (sum != 0) global_check = false;
    }
}

void check_graph(int n, ll g) { // can check either rooted or strongly connected

    int edge_cnt = __builtin_popcountll(g);
    Eigen::MatrixXd B(n, edge_cnt);

    vector<pii> edges;

    int flag = 0;
    for (int j = 0; j < n * n; j++) {
        pii e = get_edge(j+1, n); // returns 1 indexed 
        if (e.first == e.second) continue;

        if (g & (1 << j)) {
            F0R(k,n) B(k,flag) = 0; // set entire column to 0            
            B(e.first-1, flag) = 1; // -1 bc e is 1 indexed
            B(e.second-1, flag) = -1;
            flag++;
        }
        else {
            edges.pb(e);
        }
    }

    Eigen::MatrixXd B_T = B.transpose();
    Eigen::MatrixXd L = B*B_T;

    check_rowsum(n, L);

    Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(L, /* computeEigenvectors = */ false);
    auto e = es.eigenvalues(); // returns a column vector of type complex<double>

    vector<double> eigenval;
    trav(x, e) eigenval.pb(x.real());
    sort(all(eigenval));

    double base_e = eigenval[1];

    if (base_e > -eps) { // lambda_1 should be positive 
        // cout << "0: ";
    }
    else {
        // cout << endl << B << endl << L << endl;
        // cout << eigenval[1] << endl;
        global_check = false;
    }

    /*
    add 1 more edge
    */

    B.conservativeResize(B.rows(), edge_cnt+1);

    FOR(i, 0, edges.size()) { 
        F0R(j,n) B(j,flag) = 0;
        pii next = edges[i];
        B(next.first-1, flag) = 1;
        B(next.second-1, flag) = -1;

        B_T = B.transpose();
        L = B*B_T;
        es.compute(L, false);
        e = es.eigenvalues();
        eigenval.clear();
        trav(x, e) eigenval.pb(x.real());
        sort(all(eigenval));

        if (eigenval[1] - base_e >= -eps) { // lambda_1 of G' > lambda_1 of G
            // cout << "0";
        }
        else {
            // cout << endl << B << endl << L << endl;
            // cout << eigenval[1] << endl;
            global_check = false; 
        }
    }

    // cout << endl;
}

void tarjan(int v) { if (!done) {
    startTime[v] = global_time;
    low[v] = startTime[v];
    inverse_startTime[global_time] = v;
    global_time++;
    //scc_cnt++;
    stk.pb(v);

    trav(x,adj[v]) {
        if (startTime[x] == -1) {
            tarjan(x);
            low[v] = min(low[v], low[x]);
        }
        low[v] = min(low[v], startTime[x]); // do not need to check if its in another scc
    }

    if (low[v] == startTime[v]) { // identify scc
        done = true;
        while (stk[stk.size()-1] != v) {
            scc_cnt++;
            stk.pop_back();
        }
        scc_cnt++; // removes v
        stk.pop_back();
    }
} }

int check_sc(int n, ll g) {
    init();
    adj.resize(n);

    // construct adj
    for (int j = 0; j < n * n; j++) {
        if (g & (1 << j)) {
            pii e = get_edge(j+1, n); 
            adj[e.first-1].pb(e.second-1);
        }
    }

    tarjan(0); // go in from any vertex
    
    if (scc_cnt != n) return 0;
    // check_graph(n, g);
    return 1;
}

int main() {
    
    vector<ll> graphs;
    graphs.pb(2); // 1->2
    graphs.pb(6); // 1->2 , 2->1 
    
    for (int n = 3; n <= 6; n++) {
        vector<ll> new_graphs;
        int strong_cnt = 0, rooted_cnt = 0;

        trav(x,graphs) {
            ll new_g_base = 0;

            // the new graph contains all old edges
            int prev_n = n-1;
            for (int j = 0; j < prev_n * prev_n; j++) {
                if (x & (1 << j)) {
                    pii e = get_edge(j+1, prev_n); // +1 bc 1 indexed 
                    int new_e = add_edge(e.first, e.second, n);
                    new_g_base += (1 << (new_e-1)); // shift is 1 less than location
                }
            }

            // new graphs 
            for (int g = 0; g < pow(2, 2*prev_n); g++) { 
                ll new_g = new_g_base;

                // add new edges
                bool out = false;
                for (int j = 0; j < prev_n; j++) { // out edge
                    if (g & (1 << j)) {
                        int new_e = add_edge(n, j+1, n);
                        new_g += (1 << (new_e-1));
                        out = true;
                    }
                }
                bool in = false;
                for (int j = 0; j < prev_n; j++) { // in edge
                    if (g & (1 << (j + prev_n))) {
                        int new_e = add_edge(j+1, n, n);
                        new_g += (1 << (new_e)-1);
                        in = true;
                    }
                }

                // process each new graph
                if (in) {
                    // check_graph(n, new_g);
                    rooted_cnt++;
                    new_graphs.pb(new_g);
                    if (out) {// strongly connected is a subset of rooted 
                        bool tmp = check_sc(n, new_g); 
                        strong_cnt += tmp;
                    }
                }
            }
        }

        graphs = new_graphs;

        cout << "for n = " << n << endl;
        cout << global_check << endl;
        cout << "rooted: " << rooted_cnt << endl;
        cout << "strong: " << strong_cnt << endl;
        // cout << "====================\n";
    }

    return 0;
}