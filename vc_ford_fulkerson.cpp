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
#define ff first
#define ss second
#define all(x) x.begin(), x.end()
#define alla(arr, sz) arr, arr + sz
 
const int dx[4] = { 1,0,-1,0 }, dy[4] = { 0,1,0,-1 }; // for every grid problem!!
const ll linf = 4*1e18;
const int inf = 1000000007;
 
namespace io {
	void setIn(string s) { freopen(s.c_str(),"r",stdin); }
	void setOut(string s) { freopen(s.c_str(),"w",stdout); }
	void setIO(string s = "") {
		ios_base::sync_with_stdio(0); cin.tie(0); // fast I/O
		if (sz(s)) { setIn(s+".in"), setOut(s+".out"); } // for USACO
	}
}
 
using namespace io;

const int mxn = 1e5+5;

int n, root; 

vector<vi> adj, changed_adj;

vi st;
void build_st() {
    st.pb(root);
    F0R(i,n) if (i != root) {
        int l = st.size();
        int ridx = rand() % l;
        adj[st[ridx]].pb(i);
        adj[i].pb(st[ridx]);
        st.pb(i);
    }
}

void add_edges() {
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
                adj[j].pb(i);
            }
        }
    }
}

// set<int> removed;
// bool visited[20];

// void dfs(int v) {
//     visited[v] = true;
//     trav(x, adj[v]) 
//     if (removed.count(x) == 0) {
//         if (!visited[x]) {
//             dfs(x);
//         }
//     }
// }

// bool connect() {

//     F0R(t,n) {
//         if (removed.count(t) == 0 && !visited[t]) return false;
//     }
//     return true;
// }

// bool test(int k) {

//     for (int i = 0; i <= (1<<n)-1; i++) 
//     if (__builtin_popcount(i) == k) {

//         removed.clear();
//         memset(visited, false, sizeof visited);

//         for(int j = 0; j <= n-1; j++) {
//             if ((1<<j) & i) removed.insert(j);
//         }

//         F0R(j,n) {
//             if (removed.count(j) == 0) {
//                 dfs(j);
//                 if (!connect()) return false;
//                 break;
//             }
//         }
//     }

//     return true; 
// }

vector<vi> capacity;

int bfs(int s, int t, vector<int>& parent) {
    fill(parent.begin(), parent.end(), -1);
    parent[s] = -2;
    queue<pair<int, int>> q;
    q.push({s, inf});

    while (!q.empty()) {
        int cur = q.front().first;
        int flow = q.front().second;
        q.pop();

        for (int next : changed_adj[cur]) {
            if (parent[next] == -1 && capacity[cur][next]) {
                parent[next] = cur;
                int new_flow = min(flow, capacity[cur][next]);
                if (next == t)
                    return new_flow;
                q.push({next, new_flow});
            }
        }
    }

    return 0;
}

int maxflow(int s, int t) {
    int flow = 0;
    vector<int> parent((n-2)*2+2);
    int new_flow;

    while (new_flow = bfs(s, t, parent)) {
        flow += new_flow;
        int cur = t;
        while (cur != s) {
            int prev = parent[cur];
            capacity[prev][cur] -= new_flow;
            capacity[cur][prev] += new_flow;
            cur = prev;
        }
    }

    return flow;
}

#define DEBUG 0

int main() {

    srand (time(NULL));

    setOut("data.txt");

    FOR(num_nodes, 5, 100) F0R(iteration, 50) {

    #if !DEBUG

    adj.clear();
    st.clear();

    n = num_nodes;

    adj.resize(n);
    root = 0;
    build_st();

    add_edges();

    #else

    setIn("test.in");
    int edge; 
    cin >> n >> edge;
    adj.resize(n);

    F0R(i,edge) {
        int a, b;
        cin >> a >> b;
        adj[a].pb(b);
        adj[b].pb(a);
    }

    #endif 

    // F0R(i,n) {
    //     cout << i << ": "; 
    //     trav(x,adj[i]) {
    //         cout << x << ' ';
    //     }
    //     cout << '\n';
    // }

    Eigen::MatrixXd L(n,n);
    F0R(i,n) {
        F0R(j,n) L(i,j) = 0;
        L(i,i) = adj[i].size();
        trav(x,adj[i]) {
            L(i,x) = -1;
        }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(L, /* computeEigenvectors = */ false);
    auto e = es.eigenvalues(); // returns a column vector of type complex<double>

    vector<double> eigenval;
    trav(x, e) eigenval.pb(x.real());
    sort(all(eigenval));
    int l2 = eigenval[1];

    map<pii, bool> is_edge_used;

    vector<pii> edges; 
    bool mat[n][n]; memset(mat, 0, sizeof mat);
    F0R(i,n) {
        trav(x, adj[i]) {
            if (!is_edge_used[mp(i,x)]) {
                edges.pb(mp(i,x));
                is_edge_used[mp(x,i)] = true;
            }
            mat[i][x] = 1;
        }
    }

    int ans = inf;

    F0R(s,n) FOR(t,s+1,n) if (s != t && !mat[s][t]) {

        vector<pii> m(n); // map old index to new index (after splitting vertices)
        int cnt = 0;
        F0R(i,n) {
            if (s != i && t != i) {
                m[i] = mp(cnt, cnt+1);
                cnt += 2;
            }
            else {
                m[i] = mp(cnt, cnt);
                cnt++;
            }
        }

        int new_s = m[s].ff;
        int new_t = m[t].ss;

        vector<vi> adj2((n-2)*2+2);
        
        trav(x, edges) {
            pii a = m[x.ff];
            pii b = m[x.ss];

            if (a.ff == a.ss && b.ff == b.ss) {
                // ignore?  
            }

            else if (a.ff == a.ss && a.ff == new_s) {
                adj2[a.ff].pb(b.ff);
            }
            else if (a.ff == a.ss && a.ff == new_t) {
                adj2[b.ss].pb(a.ff);
            }

            else if (b.ff == b.ss && b.ff == new_s) {
                adj2[b.ff].pb(a.ff);
            }
            else if (b.ff == b.ss && b.ff == new_t) {
                adj2[a.ss].pb(b.ff);
            }

            else { // neither a nor b are source nor sink
                adj2[b.ss].pb(a.ff);
                adj2[a.ss].pb(b.ff);
            }
        }

        F0R(i,n) if (i != s && i != t) {
            pii a = m[i];
            adj2[a.ff].pb(a.ss);
        }

        capacity.clear();
        capacity.resize((n-2)*2+2);
        F0R(i,(n-2)*2+2) {
            capacity[i].resize((n-2)*2+2);
            fill(all(capacity[i]), 0);
            trav(x, adj2[i]) {
                
                if (i == new_s || new_t == x) {
                    capacity[i][x] = inf;
                }
                else {
                    capacity[i][x] = 1;
                }
            }
        }

        changed_adj.clear();
        changed_adj = adj2; 

        int res = maxflow(new_s, new_t);
        //cout << s << " " << t << " " << res << endl;
        ans = min(ans, res);
        
    }

    if (ans == inf) { // full graph
        ans = n-1;
    }

    // cout << ans << endl; 
    // cout << test(ans-1) << endl;
    // cout << test(ans) << endl;

    cout << l2 << ' ' << ans << endl; 
    
    }

    return 0;
}