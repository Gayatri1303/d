# BE_DAA_2019-PAT-SPPU


https://onlinegdb.com/uSrxQfi4R
DAA 2
#include <iostream>
#include <vector>
#include <map>
#include <ctime>
#include <algorithm>
#include <cmath>
using namespace std;

class Node {
public:
    char symbol;
    int freq;
    Node* left;
    Node* right;

    Node(char symbol, int freq, Node* left = nullptr, Node* right = nullptr)
        : symbol(symbol), freq(freq), left(left), right(right) {}
};

// Helper function to recursively calculate Huffman codes
void calculateHuffmanCodes(Node* node, const string& code, map<char, string>& huffmanCodes) {
    if (node) {
        // If it's a leaf node, add code to the map
        if (!node->left && !node->right) {
            huffmanCodes[node->symbol] = code;
        }
        // Recur for left and right children
        calculateHuffmanCodes(node->left, code + "0", huffmanCodes);
        calculateHuffmanCodes(node->right, code + "1", huffmanCodes);
    }
}

int main() {
    int n;
    cout << "Enter the number of characters: ";
    cin >> n;

    vector<char> chars(n);
    vector<int> freq(n);
    cout << "Enter characters and their frequencies:\n";
    for (int i = 0; i < n; ++i) {
        cout << "Character " << i + 1 << ": ";
        cin >> chars[i];
        cout << "Frequency of " << chars[i] << ": ";
        cin >> freq[i];
    }

    // Start measuring time for Huffman tree construction
    clock_t start_time = clock();

    // Initialize nodes list with individual character nodes
    vector<Node*> nodes;
    for (int i = 0; i < n; ++i) {
        nodes.push_back(new Node(chars[i], freq[i]));
    }

    // Build the Huffman Tree by repeatedly combining the two lowest-frequency nodes
    while (nodes.size() > 1) {
        // Sort nodes by frequency in ascending order
        sort(nodes.begin(), nodes.end(), [](Node* a, Node* b) {
            return a->freq < b->freq;
        });

        // Take the two nodes with the lowest frequencies
        Node* left = nodes[0];
        Node* right = nodes[1];

        // Create a new node combining these two
        Node* newNode = new Node('\0', left->freq + right->freq, left, right);

        // Remove the two nodes and add the new node to the list
        nodes.erase(nodes.begin());
        nodes.erase(nodes.begin());
        nodes.push_back(newNode);
    }

    // End time for Huffman tree construction
    clock_t end_time = clock();
    double tree_duration = double(end_time - start_time) / CLOCKS_PER_SEC;

    // Start measuring time for Huffman code calculation
    clock_t code_start_time = clock();

    map<char, string> huffmanCodes;
    calculateHuffmanCodes(nodes[0], "", huffmanCodes); // Root of the Huffman tree is in nodes[0]

    clock_t code_end_time = clock();
    double code_duration = double(code_end_time - code_start_time) / CLOCKS_PER_SEC;

    // Calculate estimated space in bytes for storing Huffman codes
    double spaceUsed = 0;
    for (const auto& kv : huffmanCodes) {
        int index = find(chars.begin(), chars.end(), kv.first) - chars.begin();
        spaceUsed += kv.second.length() * freq[index];
    }
    spaceUsed = ceil(spaceUsed / 8);

    // Output results
    cout << "Huffman Tree Construction Time: " << tree_duration << " seconds" << endl;
    cout << "Huffman Code Calculation Time: " << code_duration << " seconds" << endl;
    cout << "Estimated Space Used for Huffman Codes: " << spaceUsed << " bytes" << endl;

    cout << "Huffman Codes:\n";
    for (const auto& kv : huffmanCodes) {
        cout << kv.first << " -> " << kv.second << endl;
    }

    return 0;
}
DAA 1

#include <iostream>
#include <ctime> // For measuring time
using namespace std;

// Iterative Fibonacci
int fib_iterative(int n) {
    if (n <= 1) return n;
    int a = 0, b = 1, c;
    for (int i = 2; i <= n; i++) {
        c = a + b;
        a = b;
        b = c;
    }
    return b;
}

// Recursive Fibonacci
int fib_recursive(int n) {
    if (n <= 1) return n;
    return fib_recursive(n - 1) + fib_recursive(n - 2);
}

int main() {
    int n;
    cout << "Enter a number: ";
    cin >> n;

    // Measure time for iterative Fibonacci
    clock_t start_iter = clock();
    int result_iter = fib_iterative(n);
    clock_t end_iter = clock();
    double duration_iter = double(end_iter - start_iter) / CLOCKS_PER_SEC;

    // Calculate space required for iterative approach
    int space_iter = sizeof(n) + sizeof(int) * 3;  // n, a, b, c variables

    cout << "Iterative Fibonacci of " << n << " is: " << result_iter << endl;
    cout << "Time required (iterative): " << duration_iter << " seconds" << endl;
    cout << "Space required (iterative): " << space_iter << " bytes" << endl;

    // Measure time for recursive Fibonacci
    clock_t start_rec = clock();
    int result_rec = fib_recursive(n);
    clock_t end_rec = clock();
    double duration_rec = double(end_rec - start_rec) / CLOCKS_PER_SEC;

    // Calculate space required for recursive approach
    // Recursion uses a stack frame for each call, with `sizeof(int)` for each level
    int space_rec = sizeof(n) * (n > 1 ? n : 1);

    cout << "Recursive Fibonacci of " << n << " is: " << result_rec << endl;
    cout << "Time required (recursive): " << duration_rec << " seconds" << endl;
    cout << "Space required (recursive): " << space_rec << " bytes" << endl;

    return 0;
}


DAA 3

#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

// Structure to represent an item with weight and value
struct Item {
    int value, weight;

    // Constructor
    Item(int v, int w) : value(v), weight(w) {}
};

// Function to compare items based on value/weight ratio
bool compare(Item a, Item b) {
    double r1 = (double)a.value / a.weight;
    double r2 = (double)b.value / b.weight;
    return r1 > r2;  // Sort in descending order of value/weight ratio
}

// Function to solve the fractional knapsack problem
double fractionalKnapsack(int capacity, vector<Item>& items) {
    // Sort items by their value/weight ratio
    sort(items.begin(), items.end(), compare);

    double totalValue = 0.0;  // Variable to store the total value we can carry

    for (auto& item : items) {
        if (capacity >= item.weight) {
            // If the item can be completely taken, take it
            capacity -= item.weight;
            totalValue += item.value;
        } else {
            // If the item cannot be completely taken, take a fraction of it
            totalValue += item.value * ((double)capacity / item.weight);
            break;  // Knapsack is full
        }
    }

    return totalValue;
}

int main() {
    int n;  // Number of items
    int capacity;  // Knapsack capacity

    // User input
    cout << "Enter the number of items: ";
    cin >> n;

    vector<Item> items;
    int value, weight;

    // Input values and weights of items
    for (int i = 0; i < n; i++) {
        cout << "Enter value and weight of item " << i + 1 << ": ";
        cin >> value >> weight;
        items.push_back(Item(value, weight));
    }

    cout << "Enter the capacity of the knapsack: ";
    cin >> capacity;

    // Solve the fractional knapsack problem
    double maxValue = fractionalKnapsack(capacity, items);

    // Print the result
    cout << "Maximum value in knapsack = " << maxValue << endl;

    return 0;
}

DAA 5

#include <iostream>
#include <vector>
#include <ctime>
using namespace std;

bool isSafe(const vector<int>& board, int row, int col) {
    // Check previous rows for conflicts with the new queen at (row, col)
    for (int i = 0; i < row; i++) {
        int placedCol = board[i];
        // Check column conflict and both diagonals
        if (placedCol == col || abs(placedCol - col) == abs(i - row))
            return false;
    }
    return true;
}

bool solveNQueensUtil(vector<int>& board, int row, in3t n) {
    if (row == n) return true; // All queens are placed successfully

    for (int col = 0; col < n; col++) {
        if (isSafe(board, row, col)) {
            board[row] = col; // Place queen at (row, col)

            if (solveNQueensUtil(board, row + 1, n))
                return true; // Recursive call

            board[row] = -1; // Backtrack
        }
    }
    return false;
}

void printSolution(const vector<int>& board, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << (board[i] == j ? "Q " : ". ");
            cout << endl;
    }
}

int main() {
    int n;
    cout << "Enter the value of N: ";
    cin >> n;

    if (n < 1) {
        cout << "Invalid board size." << endl;
        return 1;
    }

    vector<int> board(n, -1); // Initialize board with -1 to indicate no queens placed

    clock_t start = clock();

    if (solveNQueensUtil(board, 0, n)) {
        printSolution(board, n);
    } else {
        cout << "No solution exists for N = " << n << endl;
    }

    clock_t end = clock();
    double time_taken = double(end - start) / CLOCKS_PER_SEC;
    int memory_usage = sizeof(int) * (n + n); // Rough memory estimation

    cout << "Execution time: " << time_taken << " seconds" << endl;
    cout << "Memory used: " << memory_usage << " bytes" << endl;

    return 0;
}




