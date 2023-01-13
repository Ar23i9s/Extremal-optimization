#include <iostream>
#include <random>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <numeric>
#include <igraph.h>
#include <time.h>

using namespace std;

int reverse_value(int value)
{
    if (value == 0)
        return 1;
    else
        return 0;
}

vector<vector<int>> load_matrix(const string& filepath)
{
    /// first read file, get number of lines
    ifstream data(filepath, ios::in);
    string line;
    int line_num = 0;
    while (getline(data, line))
        line_num++;

    /// read file, input data to matrix
    vector<vector<int>> matrix(line_num);
    line_num = 0;
    data = ifstream(filepath, ios::in);
    while (getline(data, line))
    {
        istringstream sin(line);
        string element;
        while (getline(sin, element, ','))
            matrix[line_num].push_back(stoi(element));
        line_num++;
    }
    return matrix;
}

void build_graph(igraph_t* graph, vector<vector<int>> sub_matrix)
{
    int sub_row_num = int(sub_matrix.size());
    int sub_col_num = int(sub_matrix[0].size());
    int node_num = sub_row_num + sub_col_num;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_empty(graph, node_num, IGRAPH_UNDIRECTED);
    igraph_vector_bool_t node_type;
    igraph_vector_bool_init(&node_type, node_num);
    for (int i = 0; i < node_num; ++i)
    {
        if (i < sub_row_num)
            VECTOR(node_type)[i] = 0;
        else
            VECTOR(node_type)[i] = 1;
        SETVAN(graph, "id", i, i);
    }

    int link_num = 0;
    for (int r = 0; r < sub_matrix.size(); ++r)
        for (int c = 0; c < sub_matrix[0].size(); ++c)
            if (sub_matrix[r][c] == 1)
                link_num++;

    int index = 0;
    igraph_vector_t edges;
    igraph_vector_init(&edges, 2 * link_num);
    for (int r = 0; r < sub_matrix.size(); ++r)
        for (int c = 0; c < sub_matrix[0].size(); ++c)
        {
            if (sub_matrix[r][c] == 1)
            {
                VECTOR(edges)[index] = r;
                VECTOR(edges)[index + 1] = c + sub_row_num;
                index += 2;
            }
        }

    igraph_create_bipartite(graph, &node_type, &edges, IGRAPH_UNDIRECTED);
    igraph_vector_bool_destroy(&node_type);
    igraph_vector_destroy(&edges);
    for (int i = 0; i < node_num; ++i)
    {
        if (i < sub_row_num)
            SETVAN(graph, "type", i, 0);
        else
            SETVAN(graph, "type", i, 1);
        SETVAN(graph, "id", i, i);
    }
}

vector<int> random_bipartition(igraph_t* graph)
{
    vector<int> row_index_list;
    vector<int> col_index_list;
    for (int i = 0; i < igraph_vcount(graph); ++i)
    {
        if (igraph_cattribute_VAN(graph, "type", i) == 0)
            row_index_list.push_back(i);
        else
            col_index_list.push_back(i);
    }

    for (int i = 0; i < col_index_list.size(); ++i)
        col_index_list[i] -= int(row_index_list.size());

    random_device rd;
    mt19937 rng(rd());
    shuffle(row_index_list.begin(), row_index_list.end(), rng);
    shuffle(col_index_list.begin(), col_index_list.end(), rng);

    int half_row_num;
    if (row_index_list.size() % 2 == 0)
        half_row_num = int(row_index_list.size() / 2);
    else
        half_row_num = int(((row_index_list.size() + 1) / 2));

    int half_col_num;
    if (col_index_list.size() % 2 == 0)
        half_col_num = int(col_index_list.size() / 2);
    else
        half_col_num = int(((col_index_list.size() + 1) / 2));

    vector<int> bipartition_row_list(row_index_list.size(), -1);
    vector<int> bipartition_col_list(col_index_list.size(), -1);
    vector<int> bipartition_list;

    for (int i = 0; i < row_index_list.size(); ++i)
    {
        if (i < half_row_num)
            bipartition_row_list[row_index_list[i]] = 0;
        else
            bipartition_row_list[row_index_list[i]] = 1;
    }

    for (int i = 0; i < col_index_list.size(); ++i)
    {
        if (i < half_col_num)
            bipartition_col_list[col_index_list[i]] = 0;
        else
            bipartition_col_list[col_index_list[i]] = 1;
    }

    bipartition_list.insert(bipartition_list.end(), bipartition_row_list.begin(), bipartition_row_list.end());
    bipartition_list.insert(bipartition_list.end(), bipartition_col_list.begin(), bipartition_col_list.end());

    return bipartition_list;
}

void adjust_bipartition(igraph_t* graph, vector<int>& bipartition)
{
    /// adjust random bipartition, change the partition of isolated nodes
    int stop_flag = 0;
    while (stop_flag == 0)
    {
        int sub_num1 = 0, sub_num2 = 0;
        for (int i : bipartition)
        {
            if (i == 0)
                sub_num1++;
            else
                sub_num2++;
        }

        igraph_vector_t index_list1;
        igraph_vector_t index_list2;
        igraph_vector_init(&index_list1, sub_num1);
        igraph_vector_init(&index_list2, sub_num2);

        sub_num1 = 0;
        sub_num2 = 0;
        for (int i = 0; i < bipartition.size(); ++i)
        {
            if (bipartition[i] == 0)
            {
                VECTOR(index_list1)[sub_num1] = i;
                sub_num1++;
            }
            else
            {
                VECTOR(index_list2)[sub_num2] = i;
                sub_num2++;
            }
        }

        igraph_t sub_graph1;
        igraph_t sub_graph2;
        igraph_vs_t vids1;
        igraph_vs_t vids2;
        igraph_vs_vector(&vids1, &index_list1);
        igraph_vs_vector(&vids2, &index_list2);
        igraph_induced_subgraph(graph, &sub_graph1, vids1, IGRAPH_SUBGRAPH_COPY_AND_DELETE);
        igraph_induced_subgraph(graph, &sub_graph2, vids2, IGRAPH_SUBGRAPH_COPY_AND_DELETE);

        stop_flag = 1;
        if (igraph_vcount(&sub_graph1) > 0)
        {
            igraph_vector_t vertex_community;
            igraph_vector_t community_size;
            igraph_vector_init(&vertex_community, 0);
            igraph_vector_init(&community_size, 0);
            igraph_integer_t cluster_num;

            igraph_clusters(&sub_graph1, &vertex_community, &community_size, &cluster_num, IGRAPH_STRONG);

            for (int i = 0; i < cluster_num; ++i)
                if (VECTOR(community_size)[i] == 1)
                {
                    long node_index;
                    igraph_vector_search(&vertex_community, 0, i, &node_index);
                    int raw_index = int(igraph_cattribute_VAN(&sub_graph1, "id", node_index));
                    bipartition[raw_index] = reverse_value(bipartition[raw_index]);
                    stop_flag = 0;
                    break;
                }
            igraph_vector_destroy(&vertex_community);
            igraph_vector_destroy(&community_size);
        }

        if (stop_flag == 1)
            if (igraph_vcount(&sub_graph2) > 0)
            {
                igraph_vector_t vertex_community;
                igraph_vector_t community_size;
                igraph_vector_init(&vertex_community, 0);
                igraph_vector_init(&community_size, 0);
                igraph_integer_t cluster_num;

                igraph_clusters(&sub_graph2, &vertex_community, &community_size, &cluster_num, IGRAPH_STRONG);

                for (int i = 0; i < cluster_num; ++i)
                    if (VECTOR(community_size)[i] == 1)
                    {
                        long node_index;
                        igraph_vector_search(&vertex_community, 0, i, &node_index);
                        int raw_index = int(igraph_cattribute_VAN(&sub_graph2, "id", node_index));
                        bipartition[raw_index] = reverse_value(bipartition[raw_index]);
                        stop_flag = 0;
                        break;
                    }
                igraph_vector_destroy(&vertex_community);
                igraph_vector_destroy(&community_size);
            }
        igraph_vector_destroy(&index_list1);
        igraph_vector_destroy(&index_list2);
        igraph_vs_destroy(&vids1);
        igraph_vs_destroy(&vids2);
        igraph_destroy(&sub_graph1);
        igraph_destroy(&sub_graph2);
    }
}

vector<int> find_community(igraph_t* graph, vector<int> bipartition)
{
    int sub_num1 = 0, sub_num2 = 0;
    for (int i : bipartition)
    {
        if (i == 0)
            sub_num1++;
        else
            sub_num2++;
    }

    igraph_vector_t index_list1;
    igraph_vector_t index_list2;
    igraph_vector_init(&index_list1, sub_num1);
    igraph_vector_init(&index_list2, sub_num2);

    sub_num1 = 0;
    sub_num2 = 0;
    for (int i = 0; i < bipartition.size(); ++i)
    {
        if (bipartition[i] == 0)
        {
            VECTOR(index_list1)[sub_num1] = i;
            sub_num1++;
        }
        else
        {
            VECTOR(index_list2)[sub_num2] = i;
            sub_num2++;
        }
    }

    igraph_t sub_graph1;
    igraph_t sub_graph2;
    igraph_vs_t vids1;
    igraph_vs_t vids2;
    igraph_vs_vector(&vids1, &index_list1);
    igraph_vs_vector(&vids2, &index_list2);
    igraph_induced_subgraph(graph, &sub_graph1, vids1, IGRAPH_SUBGRAPH_COPY_AND_DELETE);
    igraph_induced_subgraph(graph, &sub_graph2, vids2, IGRAPH_SUBGRAPH_COPY_AND_DELETE);

    /// give node the label of community
    vector<int> community_list(igraph_vcount(graph), -1);
    igraph_integer_t sub_cluster_num1 = 0;

    if (igraph_vcount(&sub_graph1) > 0)
    {
        igraph_vector_t vertex_community;
        igraph_vector_init(&vertex_community, 0);
        igraph_clusters(&sub_graph1, &vertex_community, NULL, &sub_cluster_num1, IGRAPH_STRONG);
        for (int i = 0; i < igraph_vcount(&sub_graph1); ++i)
        {
            int raw_index = int(igraph_cattribute_VAN(&sub_graph1, "id", i));
            community_list[raw_index] = int(VECTOR(vertex_community)[i]);
        }
        igraph_vector_destroy(&vertex_community);
    }

    if (igraph_vcount(&sub_graph2) > 0)
    {
        igraph_vector_t vertex_community;
        igraph_vector_init(&vertex_community, 0);
        igraph_clusters(&sub_graph2, &vertex_community, NULL, NULL, IGRAPH_STRONG);
        for (int i = 0; i < igraph_vcount(&sub_graph2); ++i)
        {
            int raw_index = int(igraph_cattribute_VAN(&sub_graph2, "id", i));
            community_list[raw_index] = int(VECTOR(vertex_community)[i] + sub_cluster_num1);
        }
        igraph_vector_destroy(&vertex_community);
    }

    igraph_vector_destroy(&index_list1);
    igraph_vector_destroy(&index_list2);
    igraph_vs_destroy(&vids1);
    igraph_vs_destroy(&vids2);
    igraph_destroy(&sub_graph1);
    igraph_destroy(&sub_graph2);

    return community_list;
}

pair<double, vector<double>> calculate_lambda_mod(vector<int> partition, vector<vector<int>> adj_matrix)
{
    int row_num = int(adj_matrix.size());
    int col_num = int(adj_matrix[0].size());
    vector<double> row_degree(row_num, 0.);
    vector<double> col_degree(col_num, 0.);

    for (int r = 0; r < row_num; ++r)
        for (int c = 0; c < col_num; ++c)
            if (adj_matrix[r][c] == 1)
            {
                row_degree[r]++;
                col_degree[c]++;
            }

    double link_num = (double)accumulate(row_degree.begin(), row_degree.end(), 0.);
    double modularity = 0.;
    vector<double> lambda_list;

    for (int i = 0; i < row_num; ++i)
    {
        double q_i = 0;
        for (int j = 0; j < col_num; ++j)
            if (partition[i] == partition[j + row_num])
                q_i += (double)adj_matrix[i][j] - row_degree[i] * col_degree[j] / link_num;
        modularity += q_i;
        lambda_list.push_back(q_i / row_degree[i]);
    }
    modularity = modularity / link_num;

    for (int j = 0; j < col_num; ++j)
    {
        double q_j = 0;
        for (int i = 0; i < row_num; ++i)
            if (partition[i] == partition[j + row_num])
                q_j += (double)adj_matrix[i][j] - row_degree[i] * col_degree[j] / link_num;
        lambda_list.push_back(q_j / col_degree[j]);
    }

    return make_pair(modularity, lambda_list);
}

pair<double, vector<double>> calculate_lambda_ibn(vector<int> partition, vector<vector<int>> adj_matrix)
{
    int row_num = int(adj_matrix.size());
    int col_num = int(adj_matrix[0].size());
    vector<double> row_degree(row_num, 0.);
    vector<double> col_degree(col_num, 0.);
    for (int r = 0; r < row_num; ++r)
        for (int c = 0; c < col_num; ++c)
            if (adj_matrix[r][c] == 1)
            {
                row_degree[r]++;
                col_degree[c]++;
            }

    map<int, vector<int>> row_index_block;
    map<int, vector<int>> col_index_block;
    for (int i = 0; i < partition.size(); ++i)
    {
        if (i < row_num)
            row_index_block[partition[i]].push_back(i);
        else
            col_index_block[partition[i]].push_back(i - row_num);
    }

    vector<int> block_index;
    for (auto index = row_index_block.begin(); index != row_index_block.end(); ++index)
        block_index.push_back(index->first);

    vector<double> lambda_list(row_num + col_num);
    for (int index : block_index)
    {
        vector<int> row_index = row_index_block[index];
        vector<int> col_index = col_index_block[index];

        /// computing row pairs overlap
        for (int j : row_index)
        {
            double lambda_value = 0.;
            double normalization = row_degree[j] * ((float)row_index.size() - 1.);
            if (normalization > 0)
            {
                for (int i : row_index)
                    if (i != j)
                    {
                        double null_model = (row_degree[i] * row_degree[j]) / (double)col_num;
                        double overlap = 0.;
                        for (int c : col_index)
                            if (adj_matrix[i][c] * adj_matrix[j][c] == 1)
                                overlap++;
                        if (row_degree[i] > row_degree[j])
                            lambda_value += overlap - null_model;
                        else if (row_degree[i] == row_degree[j])
                            lambda_value += (overlap - null_model) / 2.;
                    }
                lambda_list[j] = lambda_value / normalization;
            }
            else
                lambda_list[j] = lambda_value;
        }

        /// computing col pairs overlap
        for (int b : col_index)
        {
            double lambda_value = 0;
            double normalization = col_degree[b] * ((float)col_index.size() - 1.);
            if (normalization > 0)
            {
                for (int a : col_index)
                    if (a != b)
                    {
                        double null_model = (col_degree[a] * col_degree[b]) / (double)row_num;
                        double overlap = 0.;
                        for (int r : row_index)
                            if (adj_matrix[r][a] * adj_matrix[r][b] == 1)
                                overlap++;
                        if (col_degree[a] > col_degree[b])
                            lambda_value += overlap - null_model;
                        else if (col_degree[a] == col_degree[b])
                            lambda_value += (overlap - null_model) / 2.;
                    }
                lambda_list[b + row_num] = lambda_value / normalization;
            }
            else
                lambda_list[b + row_num] = lambda_value;
        }
    }

    double ibn = 0.;
    double lambda_sum = (double)accumulate(lambda_list.begin(), lambda_list.end(), 0.);
    ibn = (2. * lambda_sum) / (double)lambda_list.size();

    return make_pair(ibn, lambda_list);
}

pair<double, vector<double>> calculate_lambda(int objective, vector<int> partition, vector<vector<int>> adj_matrix)
{
    pair<double, vector<double>> result;
    if (objective == 0)
        result = calculate_lambda_mod(partition, adj_matrix);
    else
        result = calculate_lambda_ibn(partition, adj_matrix);
    return result;
}

double calculate_global_ibn(vector<vector<int>> adj_matrix)
{
    int row_num = int(adj_matrix.size());
    int col_num = int(adj_matrix[0].size());
    vector<double> row_degree(row_num, 0.);
    vector<double> col_degree(col_num, 0.);
    for (int r = 0; r < row_num; ++r)
        for (int c = 0; c < col_num; ++c)
            if (adj_matrix[r][c] == 1)
            {
                row_degree[r]++;
                col_degree[c]++;
            }

    vector<double> lambda_list;
    /// computing row pairs overlap
    for (int j = 0; j < row_num; ++j)
    {
        double lambda_value = 0.;
        double normalization = row_degree[j] * ((float)row_num - 1.);
        if (normalization > 0)
        {
            for (int i = 0; i < row_num; ++i)
                if (i != j)
                {
                    double null_model = (row_degree[i] * row_degree[j]) / (double)col_num;
                    double overlap = 0.;
                    for (int c = 0; c < col_num; ++c)
                        if (adj_matrix[i][c] * adj_matrix[j][c] == 1)
                            overlap++;
                    if (row_degree[i] > row_degree[j])
                        lambda_value += overlap - null_model;
                    else if (row_degree[i] == row_degree[j])
                        lambda_value += (overlap - null_model) / 2.;
                }
            lambda_list.push_back(lambda_value / normalization);
        }
        else
            lambda_list.push_back(lambda_value);
    }

    /// computing col pairs overlap
    for (int b = 0; b < col_num; ++b)
    {
        double lambda_value = 0;
        double normalization = col_degree[b] * ((float)col_num - 1.);
        if (normalization > 0)
        {
            for (int a = 0; a < col_num; ++a)
                if (a != b)
                {
                    double null_model = (col_degree[a] * col_degree[b]) / (double)row_num;
                    double overlap = 0.;
                    for (int r = 0; r < row_num; ++r)
                        if (adj_matrix[r][a] * adj_matrix[r][b] == 1)
                            overlap++;
                    if (col_degree[a] > col_degree[b])
                        lambda_value += overlap - null_model;
                    else if (col_degree[a] == col_degree[b])
                        lambda_value += (overlap - null_model) / 2.;
                }
            lambda_list.push_back(lambda_value / normalization);
        }
        else
            lambda_list.push_back(lambda_value);
    }

    double ibn = 0.;
    double lambda_sum = (double)accumulate(lambda_list.begin(), lambda_list.end(), 0.);
    ibn = (2. * lambda_sum) / (double)lambda_list.size();

    return ibn;
}

vector<int> check_movable_index(vector<int> partition, vector<vector<int>> input_matrix)
{
    /// check nodes that have neighbors in the other partition
    vector<int> movable_index;
    int row_num = int(input_matrix.size());
    int col_num = int(input_matrix[0].size());

    for (int i = 0; i < row_num; ++i)
        for (int j = 0; j < col_num; ++j)
            if (input_matrix[i][j] == 1 && partition[i] != partition[j + row_num])
            {
                movable_index.push_back(i);
                break;
            }

    for (int j = 0; j < col_num; ++j)
        for (int i = 0; i < row_num; ++i)
            if (input_matrix[i][j] == 1 && partition[i] != partition[j + row_num])
            {
                movable_index.push_back(j + row_num);
                break;
            }
    return movable_index;
}

vector<int> check_isolated_neighbors(int index, vector<int> partition, vector<vector<int>> input_matrix)
{
    /// check whether the neighbor of a node is isolated
    vector<int> isolated_neighbors;
    int row_num = int(input_matrix.size());
    int col_num = int(input_matrix[0].size());

    if (index < row_num)
    {
        for (int neighbor_index = 0; neighbor_index < col_num; ++neighbor_index)
            if (input_matrix[index][neighbor_index] == 1 && partition[index] == partition[neighbor_index + row_num])
            {
                int isolated_flag = 1;
                for (int neighbor_neighbor_index = 0; neighbor_neighbor_index < row_num; ++neighbor_neighbor_index)
                    if (input_matrix[neighbor_neighbor_index][neighbor_index] == 1 && neighbor_neighbor_index != index
                        && partition[neighbor_neighbor_index] == partition[neighbor_index + row_num])
                    {
                        isolated_flag = 0;
                        break;
                    }
                if (isolated_flag == 1)
                    isolated_neighbors.push_back(neighbor_index + row_num);
            }
    }

    else
    {
        for (int neighbor_index = 0; neighbor_index < row_num; ++neighbor_index)
            if (input_matrix[neighbor_index][index - row_num] == 1 && partition[index] == partition[neighbor_index])
            {
                int isolated_flag = 1;
                for (int neighbor_neighbor_index = 0; neighbor_neighbor_index < col_num; ++neighbor_neighbor_index)
                    if (input_matrix[neighbor_index][neighbor_neighbor_index] == 1 && neighbor_neighbor_index + row_num != index
                        && partition[neighbor_neighbor_index + row_num] == partition[neighbor_index])
                    {
                        isolated_flag = 0;
                        break;
                    }
                if (isolated_flag == 1)
                    isolated_neighbors.push_back(neighbor_index);
            }
    }

    return isolated_neighbors;
}

void proceed_move_node_ibn(int node_index, vector<double>& lambda_list, vector<int>& bipartition, igraph_t* graph,
    vector<int>& partition, vector<vector<int>> input_matrix, vector<double> row_degree, vector<double> col_degree)
{
    int row_num = int(input_matrix.size());
    int col_num = int(input_matrix[0].size());

    for (int step = 1; step < 3; ++step)
    {
        /// step 1 : move out node, modify lambda of node in old block
        /// step 2 : move in node, modify lambda of node in new block

        vector<int> row_index_list;
        vector<int> col_index_list;
        for (int i = 0; i < partition.size(); ++i)
        {
            if (partition[i] == partition[node_index])
            {
                if (i < row_num)
                    row_index_list.push_back(i);
                else
                    col_index_list.push_back(i - row_num);
            }
        }

        if (node_index < row_num)
        {
            double old_normalization;
            double new_normalization;
            if (step == 1)
            {
                old_normalization = (float)row_index_list.size() - 1.;
                new_normalization = (float)row_index_list.size() - 2.;
                lambda_list[node_index] = 0.;
            }
            else
            {
                old_normalization = (float)row_index_list.size() - 2.;
                new_normalization = (float)row_index_list.size() - 1.;

                if (old_normalization > 0)
                {
                    for (int row_index : row_index_list)
                        if (row_index != node_index && row_degree[node_index] <= row_degree[row_index])
                        {
                            double null_model = (row_degree[node_index] * row_degree[row_index]) / (double)col_num;
                            double overlap = 0.;
                            for (int col_index : col_index_list)
                                if (input_matrix[node_index][col_index] * input_matrix[row_index][col_index] == 1)
                                    overlap++;

                            if (row_degree[node_index] < row_degree[row_index])
                                lambda_list[node_index] += overlap - null_model;

                            else if (row_degree[node_index] == row_degree[row_index])
                                lambda_list[node_index] += (overlap - null_model) / 2;

                            lambda_list[node_index] = lambda_list[node_index] / (new_normalization * row_degree[node_index]);
                        }
                }
            }

            for (int row_index : row_index_list)
            {
                if ((step == 1 && new_normalization > 0) || (step == 2 && old_normalization > 0))
                {
                    if (row_index != node_index && row_degree[row_index] <= row_degree[node_index])
                    {
                        lambda_list[row_index] = lambda_list[row_index] * old_normalization;
                        double null_model = (row_degree[row_index] * row_degree[node_index]) / (double)col_num;
                        double overlap = 0.;
                        double lambda_value = 0.;
                        for (int col_index : col_index_list)
                            if (input_matrix[row_index][col_index] * input_matrix[node_index][col_index] == 1)
                                overlap++;

                        if (row_degree[row_index] < row_degree[node_index])
                            lambda_value = (overlap - null_model) / row_degree[row_index];
                        else if (row_degree[row_index] == row_degree[node_index])
                            lambda_value = (overlap - null_model) / (2. * row_degree[row_index]);

                        if (step == 1)
                            lambda_list[row_index] = (lambda_list[row_index] - lambda_value) / new_normalization;
                        else
                            lambda_list[row_index] = (lambda_list[row_index] + lambda_value) / new_normalization;
                    }
                    else
                        lambda_list[row_index] = lambda_list[row_index] * old_normalization / new_normalization;
                }
                else
                    lambda_list[row_index] = 0;
            }

            for (int col_index1 : col_index_list)
            {
                double normalization = col_degree[col_index1] * ((float)col_index_list.size() - 1.);
                if (normalization > 0)
                {
                    for (int col_index2 : col_index_list)
                    {
                        if (col_degree[col_index1] < col_degree[col_index2])
                        {
                            if (input_matrix[node_index][col_index1] * input_matrix[node_index][col_index2] == 1)
                            {
                                if (step == 1)
                                    lambda_list[col_index1 + row_num] -= 1 / normalization;
                                else
                                    lambda_list[col_index1 + row_num] += 1 / normalization;
                            }
                        }
                        else if (col_index1 != col_index2 && col_degree[col_index1] == col_degree[col_index2])
                        {
                            if (input_matrix[node_index][col_index1] * input_matrix[node_index][col_index2] == 1)
                            {
                                if (step == 1)
                                    lambda_list[col_index1 + row_num] -= 1 / (2. * normalization);
                                else
                                    lambda_list[col_index1 + row_num] += 1 / (2. * normalization);
                            }
                        }
                    }
                }
                else
                    lambda_list[col_index1 + row_num] = 0.;
            }
        }

        else
        {
            double old_normalization;
            double new_normalization;
            if (step == 1)
            {
                lambda_list[node_index] = 0.;
                old_normalization = (float)col_index_list.size() - 1.;
                new_normalization = (float)col_index_list.size() - 2.;
            }
            else
            {
                old_normalization = (float)col_index_list.size() - 2.;
                new_normalization = (float)col_index_list.size() - 1.;

                if (old_normalization > 0)
                {
                    for (int col_index : col_index_list)
                        if (col_index != node_index && col_degree[node_index - row_num] <= col_degree[col_index])
                        {
                            double null_model = (col_degree[node_index - row_num] * col_degree[col_index]) / (double)row_num;
                            double overlap = 0.;
                            for (int row_index : row_index_list)
                                if (input_matrix[row_index][node_index - row_num] * input_matrix[row_index][col_index] == 1)
                                    overlap++;

                            if (col_degree[node_index - row_num] < col_degree[col_index])
                                lambda_list[node_index] += overlap - null_model;

                            else if (col_degree[node_index - row_num] == col_degree[col_index])
                                lambda_list[node_index] += (overlap - null_model) / 2.;

                            lambda_list[node_index] = lambda_list[node_index] / (new_normalization * col_degree[node_index - row_num]);
                        }
                }
            }

            for (int col_index : col_index_list)
            {
                if ((step == 1 && new_normalization > 0) || (step == 2 && old_normalization > 0))
                {
                    if (col_index != node_index - row_num && col_degree[col_index] <= col_degree[node_index - row_num])
                    {
                        lambda_list[col_index + row_num] = lambda_list[col_index + row_num] * old_normalization;
                        double null_model = (col_degree[col_index] * col_degree[node_index - row_num]) / (double)row_num;
                        double overlap = 0.;
                        double lambda_value = 0.;
                        for (int row_index : row_index_list)
                            if (input_matrix[row_index][col_index] * input_matrix[row_index][node_index - row_num] == 1)
                                overlap++;

                        if (col_degree[col_index] < col_degree[node_index - row_num])
                            lambda_value = (overlap - null_model) / col_degree[col_index];
                        else if (col_degree[col_index] == col_degree[node_index - row_num])
                            lambda_value = (overlap - null_model) / (2. * col_degree[col_index]);

                        if (step == 1)
                            lambda_list[col_index + row_num] = (lambda_list[col_index + row_num] - lambda_value) / new_normalization;
                        else
                            lambda_list[col_index + row_num] = (lambda_list[col_index + row_num] + lambda_value) / new_normalization;
                    }

                    else
                        lambda_list[col_index + row_num] = lambda_list[col_index + row_num] * old_normalization / new_normalization;
                }
                else
                    lambda_list[col_index + row_num] = 0;
            }

            for (int row_index1 : row_index_list)
            {
                double normalization = row_degree[row_index1] * ((float)row_index_list.size() - 1.);
                if (normalization > 0)
                {
                    for (int row_index2 : row_index_list)
                    {
                        if (row_degree[row_index1] < row_degree[row_index2])
                        {
                            if (input_matrix[row_index1][node_index - row_num] * input_matrix[row_index2][node_index - row_num] == 1)
                            {
                                if (step == 1)
                                    lambda_list[row_index1] -= 1 / normalization;
                                else
                                    lambda_list[row_index1] += 1 / normalization;
                            }
                        }
                        else if (row_index1 != row_index2 && row_degree[row_index1] == row_degree[row_index2])
                        {
                            if (input_matrix[row_index1][node_index - row_num] * input_matrix[row_index2][node_index - row_num] == 1)
                            {
                                if (step == 1)
                                    lambda_list[row_index1] -= 1 / (2. * normalization);
                                else
                                    lambda_list[row_index1] += 1 / (2. * normalization);
                            }
                        }
                    }
                }
                else
                    lambda_list[row_index1] = 0.;
            }
        }

        if (step == 1)
        {
            bipartition[node_index] = reverse_value(bipartition[node_index]);
            partition = find_community(graph, bipartition);
        }
    }
}

void update_community(vector<int>& community_list)
{
    ///:return Sequential community index from 0
    int community_index = 0;
    while (community_index < *max_element(community_list.begin(), community_list.end()))
    {
        if (find(community_list.begin(), community_list.end(), community_index) != community_list.end())
            community_index++;
        else
        {
            for (int i = 0; i < community_list.size(); ++i)
                if (community_list[i] > community_index)
                    community_list[i]--;
        }
    }
}

void print_community(int row_num, int col_num, vector<int> community_list)
{
    cout << "community: ";
    for (int r : community_list)
        cout << r << ", ";
    cout << endl;

    vector<int> row_community(community_list.begin(), community_list.begin() + row_num);
    sort(row_community.begin(), row_community.end());
    vector<int>::iterator pos = unique(row_community.begin(), row_community.end());
    row_community.erase(pos, row_community.end());
    cout << "row: ";
    for (int r : row_community)
        cout << r << " ";
    cout << endl;

    vector<int> col_community(community_list.begin() + row_num, community_list.end());
    sort(col_community.begin(), col_community.end());
    pos = unique(col_community.begin(), col_community.end());
    col_community.erase(pos, col_community.end());
    cout << "col: ";
    for (int c : col_community)
        cout << c << " ";
    cout << endl;
}

vector<int> call_extremal_optimization(vector<vector<int>> sub_matrix, int objective, double alpha_ratio)
{
    int sub_row_num = int(sub_matrix.size());
    int sub_col_num = int(sub_matrix[0].size());

    vector<double> sub_row_degree(sub_row_num, 0.);
    vector<double> sub_col_degree(sub_col_num, 0.);
    for (int r = 0; r < sub_row_num; ++r)
        for (int c = 0; c < sub_col_num; ++c)
            if (sub_matrix[r][c] == 1)
            {
                sub_row_degree[r]++;
                sub_col_degree[c]++;
            }

    double metric_best = 0.;
    vector<int> sub_community_best(sub_row_num + sub_col_num, 0);

    /// build bipartite graph
    igraph_t graph;
    build_graph(&graph, sub_matrix);

    /// initial bipartition
    vector<int> bipartition_list = random_bipartition(&graph);
    adjust_bipartition(&graph, bipartition_list);
    vector<int> sub_community = find_community(&graph, bipartition_list);

    double metric_new = 0.;
    vector<double> lambda_list(sub_row_num + sub_col_num);
    pair<double, vector<double>> result = calculate_lambda(objective, sub_community, sub_matrix);
    metric_new = result.first;
    lambda_list = result.second;

    metric_best = metric_new;
    copy(sub_community.begin(), sub_community.end(), sub_community_best.begin());

    /// optimization
    int alpha = 0;
    int alpha_max = int((sub_row_num + sub_col_num) * alpha_ratio);
    while (alpha < alpha_max)
    {
        /// select and move a node to the other partition
        vector<int> movable_index = check_movable_index(sub_community, sub_matrix);
        vector<int> sub_community_copy(sub_row_num + sub_col_num, 0);
        vector<int> bipartition_list_copy(sub_row_num + sub_col_num, 0);
        vector<double> lambda_list_copy(sub_row_num + sub_col_num, 0);

        if (movable_index.size() > 0)
        {
            copy(sub_community.begin(), sub_community.end(), sub_community_copy.begin());
            copy(bipartition_list.begin(), bipartition_list.end(), bipartition_list_copy.begin());
            copy(lambda_list.begin(), lambda_list.end(), lambda_list_copy.begin());

            sort(movable_index.begin(), movable_index.end(),
                [&lambda_list](int x, int y) {return lambda_list[x] < lambda_list[y]; });

            double tau = 1. + 1. / log(sub_row_num + sub_col_num);
            vector<double> prob_list;
            for (int i = 0; i < movable_index.size(); ++i)
                prob_list.push_back(pow((i + 1), -tau));

            random_device rd;
            mt19937 gen(rd());
            discrete_distribution<int> distribution(prob_list.begin(), prob_list.end());
            int select_index = movable_index[distribution(gen)];

            vector<int> move_index = { select_index };
            vector<int> follow_index = check_isolated_neighbors(select_index, sub_community, sub_matrix);
            for (int index : follow_index)
                move_index.push_back(index);

            if (objective == 0)
            {
                for (int index : move_index)
                    bipartition_list[index] = reverse_value(bipartition_list[index]);
                sub_community = find_community(&graph, bipartition_list);
                result = calculate_lambda(objective, sub_community, sub_matrix);
                metric_new = result.first;
                lambda_list = result.second;
            }
            else
            {
                for (int index : move_index)
                    proceed_move_node_ibn(index, lambda_list, bipartition_list, &graph,
                        sub_community, sub_matrix, sub_row_degree, sub_col_degree);
                double lambda_sum = (double)accumulate(lambda_list.begin(), lambda_list.end(), 0.);
                metric_new = (2. * lambda_sum) / (double)lambda_list.size();
            }

            if (metric_new > metric_best)
            {
                copy(sub_community.begin(), sub_community.end(), sub_community_best.begin());
                metric_best = metric_new;
                alpha = 0;
            }
            else
            {
                copy(bipartition_list_copy.begin(), bipartition_list_copy.end(), bipartition_list.begin());
                copy(lambda_list_copy.begin(), lambda_list_copy.end(), lambda_list.begin());
                copy(sub_community_copy.begin(), sub_community_copy.end(), sub_community.begin());
                alpha++;
            }
        }

        else
            break;
    }

    update_community(sub_community_best);
    return sub_community_best;
}

pair<double, vector<int>> extremal_optimization(const string& filepath, int objective, double alpha_ratio, int iteration)
{
    /// objective : modularity - 0, in-block nestedness - 1
    vector<vector<int>> input_matrix = load_matrix(filepath);

    int row_num = int(input_matrix.size());
    int col_num = int(input_matrix[0].size());
    double metric_best_step = 0.;
    vector<int> community_best_step(row_num + col_num, 0);
    double metric_best_global = 0.;
    vector<int> community_best_global(row_num + col_num, 0);
    pair<double, vector<double>> result;

    for (int i = 0; i < iteration; ++i)
    {
        int community_num = *max_element(community_best_step.begin(), community_best_step.end()) + 1;
        vector<vector<int>> row_index_community(community_num);
        vector<vector<int>> col_index_community(community_num);
        vector<vector<int>> raw_index_community(community_num);

        for (int community_index = 0; community_index < community_num; ++community_index)
            for (int index = 0; index < community_best_step.size(); ++index)
                if (community_best_step[index] == community_index)
                {
                    if (index < row_num)
                        row_index_community[community_index].push_back(index);
                    else
                        col_index_community[community_index].push_back(index - row_num);
                    raw_index_community[community_index].push_back(index);
                }

        for (int community_index = 0; community_index < community_num; ++community_index)
        {
            int sub_row_num = int(row_index_community[community_index].size());
            int sub_col_num = int(col_index_community[community_index].size());
            vector<vector<int>> sub_matrix(sub_row_num, vector<int>(sub_col_num));

            if (sub_row_num > 1 && sub_col_num > 1)
            {
                for (int r = 0; r < sub_row_num; ++r)
                    for (int c = 0; c < sub_col_num; ++c)
                        sub_matrix[r][c] = input_matrix[row_index_community[community_index][r]][col_index_community[community_index][c]];

                vector<int> sub_community_best = call_extremal_optimization(sub_matrix, objective, alpha_ratio);
                int max_community_index = *max_element(community_best_step.begin(), community_best_step.end());
                for (int index = 0; index < sub_community_best.size(); ++index)
                    community_best_step[raw_index_community[community_index][index]] = sub_community_best[index] + max_community_index + 1;
                result = calculate_lambda(objective, community_best_step, input_matrix);
                metric_best_step = result.first;
                update_community(community_best_step);

                if (metric_best_step > metric_best_global)
                {
                    metric_best_global = metric_best_step;
                    copy(community_best_step.begin(), community_best_step.end(), community_best_global.begin());
                }
                else
                    copy(community_best_global.begin(), community_best_global.end(), community_best_step.begin());
            }
        }
    }

    printf("metric_best_global %f\n", metric_best_global);
    print_community(row_num, col_num, community_best_global);

    if (metric_best_global <= 0)
        return make_pair(calculate_global_ibn(input_matrix), vector<int>(row_num + col_num, 0));
    else
        return make_pair(metric_best_global, community_best_global);
}

void start_save(const string& filepath, int objective, double alpha_ratio, int iteration, int loop)
{
    string basepath = filepath.substr(0, filepath.find_last_of("\\") + 1);
    string filename = filepath.substr(filepath.find_last_of("\\") + 1,
        filepath.find_last_of(".") - filepath.find_last_of("\\") - 1);
    vector<vector<int>> matrix = load_matrix(filepath);

    ofstream data;
    if (objective == 1)
    {
        string appendix = "_" + to_string(alpha_ratio).substr(0, 3) + "_" + to_string(iteration) + "_ibn.txt";
        filename += appendix;
        data.open(basepath + filename, ios::app);
        double global_ibn = calculate_global_ibn(matrix);
        data << global_ibn << endl;
    }

    else
    {
        string appendix = "_" + to_string(alpha_ratio).substr(0, 3) + "_" + to_string(iteration) + "_mod.txt";
        filename += appendix;
        data.open(basepath + filename, ios::app);
    }

    for (int i = 0; i < loop; ++i)
    {
        clock_t start = clock();
        pair<double, vector<int>> result = extremal_optimization(filepath, objective, alpha_ratio, iteration);
        data << result.first << endl;
        for (int i = 0; i < result.second.size(); ++i)
        {
            if (i != matrix.size() - 1 && i != result.second.size() - 1)
                data << result.second[i] << ", ";
            else
                data << result.second[i] << endl;
        }

        data << "running time " << (double)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    }

    data.close();
}

int main(int argc, char* argv[])
{
    start_save(argv[1], atoi(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]));
    return 0;
}
