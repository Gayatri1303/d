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

