#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <random>
#include <unordered_set>
#include <algorithm>
using namespace std;


class MaxCliqueTabuSearch
{
public:
    static int GetRandom(int a, int b)
    {
        static mt19937 generator;
        uniform_int_distribution<int> uniform(a, b);
        return uniform(generator);
    }

    void ReadGraphFile(string filename)
    {
        ifstream fin(filename);
        string line;
        int vertices = 0, edges = 0;
        while (getline(fin, line))
        {
            if (line[0] == 'c')
            {
                continue;
            }

            stringstream line_input(line);
            char command;
            if (line[0] == 'p')
            {
                string type;
                line_input >> command >> type >> vertices >> edges;
                neighbour_sets.resize(vertices);
                qco.resize(vertices);
                index.resize(vertices, -1);
                non_neighbours.resize(vertices);
            }
            else
            {
                int start, finish;
                line_input >> command >> start >> finish;
                // Edges in DIMACS file can be repeated, but it is not a problem for our sets
                neighbour_sets[start - 1].insert(finish - 1);
                neighbour_sets[finish - 1].insert(start - 1);
            }
        }
        for (int i = 0; i < vertices; ++i)
        {
            for (int j = 0; j < vertices; ++j)
            {
                if (neighbour_sets[i].count(j) == 0 && i != j)
                    non_neighbours[i].insert(j);
            }
        }
    }

    void RunSearch(int starts, int randomization)
    {
        for (int iter = 0; iter < starts; ++iter)
        {
            ClearClique();
            for (size_t i = 0; i < neighbour_sets.size(); ++i)
            {
                qco[i] = i;
                index[i] = i;
            }
            RunInitialHeuristic(randomization);
            c_border = q_border;
            int swaps = 0;
            while (swaps < 100)
            {
                if (! Move())
                {
                    if (! Swap1To1())
                    {
                        break;
                    }
                    else
                    {
                        ++swaps;
                    }
                }
            }
            if (q_border > best_clique.size())
            {
                best_clique.clear();
                for (int i = 0; i < q_border; ++i)
                    best_clique.insert(qco[i]);
            }
        }
    }

    const unordered_set<int>& GetClique()
    {
        return best_clique;
    }

    bool Check()
    {
        for (int i : best_clique)
        {
            for (int j : best_clique)
            {
                if (i != j && neighbour_sets[i].count(j) == 0)
                {
                    cout << "Returned subgraph is not a clique\n";
                    return false;
                }
            }
        }
        return true;
    }

    void ClearClique()
    {
        best_clique.clear();
        q_border = 0;
        c_border = 0;
    }

private:
    int ComputeTightness(int vertex)
    {
        int tightness = 0;
        for (int i = 0; i < q_border; ++i)
        {
            if (neighbour_sets[qco[i]].count(vertex) == 0)
                ++tightness;
        }
        return tightness;
    }

    void SwapVertices(int vertex, int border)
    {
        int vertex_at_border = qco[border];
        swap(qco[index[vertex]], qco[border]);
        swap(index[vertex], index[vertex_at_border]);
    }

    void InsertToClique(int i)
    {
        for (int j : non_neighbours[i])
        {
            if (ComputeTightness(j) == 0)
            {
                --c_border;
                SwapVertices(j, c_border);
            }
        }
        SwapVertices(i, q_border);
        ++q_border;
    }

    void RemoveFromClique(int k)
    {
        for (int j : non_neighbours[k])
        {
            if (ComputeTightness(j) == 1)
            {
                SwapVertices(j, c_border);
                c_border++;
            }
        }
        --q_border;
        SwapVertices(k, q_border);
    }

    bool Swap1To1()
    {
        int st = GetRandom(0, q_border - 1);
        for (int counter = 0; counter < q_border; ++counter)
        {
            int vertex_index = (counter + st) % q_border;
            int vertex = qco[vertex_index];
            vector<int> L;
            for (int i : non_neighbours[vertex])
            {
                if (ComputeTightness(i) == 1)
                {
                    L.push_back(i);
                }
            }
            if (L.empty())
                continue;
            int index_in_l = GetRandom(0, L.size() - 1);
            int change = L[index_in_l];
            RemoveFromClique(vertex);
            InsertToClique(change);
            return true;
        }
        return false;
    }

    bool Move()
    {
        if (c_border == q_border)
            return false;
        int index_in_qco = GetRandom(q_border, c_border - 1);
        int vertex = qco[index_in_qco];
        InsertToClique(vertex);
        return true;
    }

    void RunInitialHeuristic(int randomization)
    {
        static mt19937 generator;
        vector<int> candidates(neighbour_sets.size());
        for (size_t i = 0; i < neighbour_sets.size(); ++i)
        {
            candidates[i] = i;
        }
        shuffle(candidates.begin(), candidates.end(), generator);
        while (! candidates.empty())
        {
            int last = candidates.size() - 1;
            int rnd = GetRandom(0, min(randomization - 1, last));
            int vertex = candidates[rnd];
            SwapVertices(vertex, q_border);
            ++q_border;
            for (int c = 0; c < candidates.size(); ++c)
            {
                int candidate = candidates[c];
                if (neighbour_sets[vertex].count(candidate) == 0)
                {
                    // Move the candidate to the end and pop it
                    swap(candidates[c], candidates[candidates.size() - 1]);
                    candidates.pop_back();
                    --c;
                }
            }
            shuffle(candidates.begin(), candidates.end(), generator);
        }
    }

private:
    vector<unordered_set<int>> neighbour_sets;
    vector<unordered_set<int>> non_neighbours;
    unordered_set<int> best_clique;
    vector<int> qco;
    vector<int> index;
    int q_border = 0;
    int c_border = 0;
};

class MaxCliqueProblem
{
public:
    // Генерирование случайных чисел на отрезке [a, b]
    static int GetRandom(int a, int b)
    {
        static mt19937 generator(time(nullptr)); // Сбросим состояние генератора (по времени)
        uniform_int_distribution<int> uniform(a, b);
        return uniform(generator);
    }

    // Чтение графа из файла формата DIMACS
    void ReadGraphFile(string filename)
    {
        ifstream fin(filename);
        string line;
        int vertices = 0, edges = 0;
        while (getline(fin, line))
        {
            if (line[0] == 'c') continue; // Пропускаем комментарии

            stringstream line_input(line);
            char command;
            if (line[0] == 'p') // Обрабатываем строку с параметрами
            {
                string type;
                line_input >> command >> type >> vertices >> edges;
                neighbour_sets.resize(vertices);
            }
            else // Добавляем ребра
            {
                int start, finish;
                line_input >> command >> start >> finish;
                neighbour_sets[start - 1].insert(finish - 1);
                neighbour_sets[finish - 1].insert(start - 1);
            }
        }
    }

    // Реализация жадного рандомизированного поиска максимальной клики
    void GreedyRandomizedMaximumClique(int randomization, int iterations)
    {
        static mt19937 generator(time(nullptr)); // Локальный экземпляр генерации
        for (int iteration = 0; iteration < iterations; ++iteration)
        {
            vector<int> current_clique; // Текущая найденная клика
            vector<int> candidates(neighbour_sets.size());
            // Заполняем кандидатов вершинами
            for (int i = 0; i < neighbour_sets.size(); ++i)
            {
                candidates[i] = i;
            }
            // Перемешиваем кандидатов
            shuffle(candidates.begin(), candidates.end(), generator);

            // Построение клики жадным образом
            while (!candidates.empty())
            {
                int last = candidates.size() - 1;
                int rnd = GetRandom(0, min(randomization - 1, last)); // Случайный выбор
                int vertex = candidates[rnd]; // Берем вершину
                current_clique.push_back(vertex);

                // Фильтруем список кандидатов, оставляя только соседей
                for (int c = 0; c < candidates.size(); ++c)
                {
                    int candidate = candidates[c];
                    if (neighbour_sets[vertex].count(candidate) == 0) // Если нет ребра
                    {
                        // Убираем кандидата (обмениваем с последним)
                        swap(candidates[c], candidates[last]);
                        candidates.pop_back();
                        --c;
                        --last;
                    }
                }
                shuffle(candidates.begin(), candidates.end(), generator); // Перемешиваем остаток
            }
            // Сохраняем лучшую клику, если текущая оказалась больше
            if (current_clique.size() > best_clique.size())
            {
                best_clique = current_clique;
            }
        }
    }

    // Возвращает найденную максимальную клику
    unordered_set<int> GetClique() const
    {
        unordered_set<int> best_clique2;
        for (const int& element : best_clique) // Если best_clique — это std::vector<int> или любой другого типа контейнер
        {
            best_clique2.insert(element);
        }
        return best_clique2; // Возвращаем копию множества
    }

    // Проверка найденной клики на корректность
    bool Check()
    {
        // Проверяем, что клика не содержит дубликатов
        if (unique(best_clique.begin(), best_clique.end()) != best_clique.end())
        {
            cout << "Duplicated vertices in the clique\n";
            return false;
        }
        // Проверяем, что каждая пара вершин в клике соединена ребром
        for (int i : best_clique)
        {
            for (int j : best_clique)
            {
                if (i != j && neighbour_sets[i].count(j) == 0)
                {
                    cout << "Returned subgraph is not a clique\n";
                    return false;
                }
            }
        }
        return true;
    }

private:
    vector<unordered_set<int>> neighbour_sets; // Список соседей для каждой вершины
    vector<int> best_clique; // Лучшая найденная клика
};


class BnBSolver
{

    static int GetRandom(int a, int b)
    {
        static mt19937 generator;
        uniform_int_distribution<int> uniform(a, b);
        return uniform(generator);
    }


public:
    void ReadGraphFile(string filename)
    {
        ifstream fin(filename);
        string line;
        file = filename;
        int vert = 0, edges = 0;
        while (getline(fin, line))
        {
            if (line[0] == 'c')
            {
                continue;
            }
            if (line[0] == 'p')
            {
                stringstream s(line);
                char c;
                string in;
                s >> c >> in >> vert >> edges;
                neighbours.resize(vert);
            }
            else
            {
                stringstream s(line);
                char c;
                int st, fn;
                s >> c >> st >> fn;
                neighbours[st - 1].insert(fn - 1);
                neighbours[fn - 1].insert(st - 1);
            }
        }
    }

    void SortVerticesByDegree() {
        vector<size_t> indices(neighbours.size());
        iota(indices.begin(), indices.end(), 0);
        sort(indices.begin(), indices.end(), [&](int u, int v) {
            return neighbours[u].size() < neighbours[v].size();
        });

        vector<unordered_set<int>> new_neighbours(indices.size());
        for (size_t i = 0; i < indices.size(); i++) {
            for (int neighbor : neighbours[indices[i]]) {
                int new_idx = find(indices.begin(), indices.end(), neighbor) - indices.begin();
                new_neighbours[i].insert(new_idx);
            }
        }
        neighbours = std::move(new_neighbours);
    }

    void PreprocessGraph() {
        int min_degree = best_clique.size();
        vector<int> filtered_vertices;
        for (size_t i = 0; i < neighbours.size(); i++) {
            if (neighbours[i].size() >= min_degree) {
                filtered_vertices.push_back(i);
            }
        }

        vector<unordered_set<int>> new_neighbours(filtered_vertices.size());
        for (size_t i = 0; i < filtered_vertices.size(); i++) {
            for (int neighbor : neighbours[filtered_vertices[i]]) {
                if (binary_search(filtered_vertices.begin(), filtered_vertices.end(), neighbor)) {
                    new_neighbours[i].insert(neighbor);
                }
            }
        }
        neighbours = std::move(new_neighbours);
    }

    void RunBnB()
    {
        MaxCliqueProblem st;
        st.ReadGraphFile(file);
        st.GreedyRandomizedMaximumClique(10000, 10000);
        best_clique = st.GetClique();

        MaxCliqueTabuSearch st2;
        st2.ReadGraphFile(file);
        st2.RunSearch(1, 10);
        if (st2.GetClique().size() > best_clique.size())
            best_clique = st2.GetClique();



        PreprocessGraph();     // Уменьшаем размер графа
        SortVerticesByDegree(); // Сортируем для лучшего ветвления

        vector<int> candidates(neighbours.size());
        for (size_t i = 0; i < neighbours.size(); ++i)
        {
            candidates[i] = i;
        }
        // static mt19937 generator;
        // shuffle(candidates.begin(), candidates.end(), generator);
        BnBWithNodeRemoval(candidates);
    }

    const unordered_set<int>& GetClique()
    {
        return best_clique;
    }

    bool Check()
    {
        for (int i : clique)
        {
            for (int j : clique)
            {
                if (i != j && neighbours[i].count(j) == 0)
                {
                    cout << "Returned subgraph is not clique\n";
                    return false;
                }
            }
        }
        return true;
    }

    void ClearClique()
    {
        best_clique.clear();
        clique.clear();
    }

private:

    // Определение списка узлов, связанных со всеми, кроме одного узла в текущей клике
    vector<int> MakeOneMissing()
    {
        vector<int> result;
        for (int i = 0; i < neighbours.size(); ++i)
        {
            if (clique.count(i)) // Узел уже в клике
                continue;

            int count = 0;
            for (int v : clique)
            {
                if (neighbours[i].count(v))
                    ++count;
            }

            // Узел подключен ко всем, кроме одного
            if (count == clique.size() - 1)
                result.push_back(i);
        }
        return result;
    }

    // Определение узла, который лучше всего удалить
    int GetNodeToDrop(const vector<int>& oneMissing) {
        int node_to_drop = -1;
        int max_degree_in_clique = 0;

        for (int node : clique) {
            int count_connections_in_missing = 0;
            for (int m : oneMissing) {
                if (neighbours[node].count(m)) {
                    ++count_connections_in_missing;
                }
            }

            if (count_connections_in_missing >= max_degree_in_clique) {
                max_degree_in_clique = count_connections_in_missing;
                node_to_drop = node;
            }
        }

        return node_to_drop; // Лучший узел для удаления
    }

    // Основной метод рекурсии для поиска клики с учетом удаления узлов
    void BnBWithNodeRemoval(vector<int>& candidates)
    {
        // Если кандидатов нет, обновляем лучшую клику
        if (candidates.empty())
        {
            if (clique.size() > best_clique.size())
                best_clique = clique;
            return;
        }

        // Улучшение: если клика и кандидаты не могут превзойти текущую лучшую
        if (clique.size() + candidates.size() <= best_clique.size())
            return;

        for (size_t c = 0; c < candidates.size(); ++c)
        {
            // Улучшение: если клика и кандидаты не могут превзойти текущую лучшую
            // if (clique.size() + candidates.size() - c + 1 <= best_clique.size())
            //     return;
            int current_vertex = candidates[c];
            vector<int> new_candidates;
            for (size_t i = c + 1; i < candidates.size(); ++i)
            {
                if (neighbours[current_vertex].count(candidates[i])) // Являются соседями
                    new_candidates.push_back(candidates[i]);
            }

            // Добавляем текущую вершину в клику
            clique.insert(current_vertex);

            // Если кандидатов нет, обновляем лучшую клику
            if (new_candidates.empty())
            {
                if (clique.size() > best_clique.size())
                    best_clique = clique;

                // Убираем текущую вершину из клики
                clique.erase(current_vertex);
                continue;
            }

            // Улучшение: если клика и кандидаты не могут превзойти текущую лучшую
            if (clique.size() + new_candidates.size() <= best_clique.size()) {
                // Убираем текущую вершину из клики
                clique.erase(current_vertex);
                continue;
            }


            sort(new_candidates.begin(), new_candidates.end(), [this](int u, int v) {
                int u_intersection = 0, v_intersection = 0;
                for (int w : clique) {
                    if (neighbours[u].count(w)) ++u_intersection;
                    if (neighbours[v].count(w)) ++v_intersection;
                }
                return (u_intersection < v_intersection) ||
                       (u_intersection == v_intersection && neighbours[u].size() < neighbours[v].size());
            });



            // Улучшение с помощью удаления узлов
            if (!clique.empty())
            {
                vector<int> oneMissing = MakeOneMissing();
                int node_to_drop = GetNodeToDrop(oneMissing);
                if (node_to_drop != -1)
                {
                    clique.erase(node_to_drop);

                    // Обновляем список кандидатов
                    vector<int> updated_candidates;
                    for (int candidate : new_candidates)
                    {
                        if (neighbours[node_to_drop].count(candidate) == 0 || clique.count(candidate) > 0)
                            updated_candidates.push_back(candidate);
                    }

                    BnBWithNodeRemoval(updated_candidates);

                    // Восстанавливаем клику
                    clique.insert(node_to_drop);
                }
            }

            // Убираем текущую вершину из клики
            clique.erase(current_vertex);
        }
    }


private:
    vector<unordered_set<int>> neighbours;
    unordered_set<int> best_clique;
    unordered_set<int> clique;
    string file;
};

int main()
{
    // ios_base::sync_with_stdio(false);
    // cin.tie(nullptr);
    vector<string> files = { //"brock200_1.clq",
        "brock200_2.clq", "brock200_3.clq", "brock200_4.clq",
        "C125.9.clq",
        "gen200_p0.9_44.clq", "gen200_p0.9_55.clq",
        "hamming8-4.clq",
       "johnson16-2-4.clq", "johnson8-2-4.clq",
        "keller4.clq",
        "MANN_a27.clq", "MANN_a9.clq",
        "p_hat1000-1.clq",  "p_hat1500-1.clq", "p_hat300-3.clq",
        "san1000.clq",
        "sanr200_0.9.clq" };
    // ofstream fout("clique_bnb.csv");
    // fout << "File; Clique; Time (sec)\n";
    for (string file : files)
    {
        cout << file;
        BnBSolver problem;
        problem.ReadGraphFile(file);
        problem.ClearClique();
        clock_t start = clock();
        problem.RunBnB();
        if (! problem.Check())
        {
            cout << "*** WARNING: incorrect clique ***\n";
            // fout << "*** WARNING: incorrect clique ***\n";
        }
        // fout << file << "; " << problem.GetClique().size() << "; " << double(clock() - start) / 1000 << '\n';
        cout << ", result - " << problem.GetClique().size() << ", time - " << double(clock() - start) / 1000 << '\n';
    }
    return 0;
}
