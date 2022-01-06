//
// Created by Ryan.Zurrin001 on 1/5/2022.
//
#include<bits/stdc++.h>

#ifndef STATBST_T_BST_H
#define STATBST_T_BST_H

using namespace std;
template<class T>
class t_treeNode {
public:
    [[maybe_unused]] T data;
    [[maybe_unused]] t_treeNode<T> *left;
    [[maybe_unused]] t_treeNode<T> *right;
    [[maybe_unused]] t_treeNode<T> *parent;
    [[maybe_unused]] explicit t_treeNode(T data) {
        this->data = data;
        left = nullptr;
        right = nullptr;
        parent = nullptr;
    }
}; // end of class t_treeNode

template <class T>
class [[maybe_unused]] t_BST {
private:
    // private instance variables
    [[maybe_unused]] t_treeNode<T> *root;
    [[maybe_unused]] int size{};
    // private insert methods
    [[maybe_unused]] void insert(t_treeNode<T> *&root_, t_treeNode<T> *&newNode);
    // private delete methods
    [[maybe_unused]] void remove(t_treeNode<T> *&root_, T data);
    [[maybe_unused]] void removeAll(t_treeNode<T> *&root_, T data);
    // private print methods
    [[maybe_unused]] void printTree(t_treeNode<T> *root_);
    [[maybe_unused]] void inOrder(t_treeNode<T> *root_);
    [[maybe_unused]] void preOrder(t_treeNode<T> *root_);
    [[maybe_unused]] void postOrder(t_treeNode<T> *node);
    [[maybe_unused]] void levelOrder(t_treeNode<T> *node);
    [[maybe_unused]] void print2DUtil(t_treeNode<T> *root_, int space);
    // private search methods
    [[maybe_unused]] t_treeNode<T> *findSuccessor(t_treeNode<T> *node);
    [[maybe_unused]] t_treeNode<T> *findPredecessor(t_treeNode<T> *node);
    [[maybe_unused]] t_treeNode<T> *find(t_treeNode<T> *node, T data);
    [[maybe_unused]] t_treeNode<T> *findMin(t_treeNode<T> *node);
    [[maybe_unused]] t_treeNode<T> *findMax(t_treeNode<T> *node);
    // private tree balance methods
    [[maybe_unused]] void rightRotate(t_treeNode<T> *&node);
    [[maybe_unused]] void leftRotate(t_treeNode<T> *&node);
    [[maybe_unused]] void leftRightRotate(t_treeNode<T> *&node);
    [[maybe_unused]] void rightLeftRotate(t_treeNode<T> *&node);
    [[maybe_unused]] void balance(t_treeNode<T> *&node);
    [[maybe_unused]] int getBalanceFactor(t_treeNode<T> *node);
    // private getter methods
    [[maybe_unused]] int getLeaves(t_treeNode<T> *node);
    [[maybe_unused]] int getHeight(t_treeNode<T> *node);
    [[maybe_unused]] int getDiameter(t_treeNode<T>* root_);
    [[maybe_unused]] int getCount(t_treeNode<T> *node);
    [[maybe_unused]] int keyCountUtil(t_treeNode<T>* root_, T key);
    // private transformation methods
    [[maybe_unused]] T replaceWithSums(t_treeNode<T> *&node);
    // private static methods
    [[maybe_unused]] double sum(t_treeNode<T> *node);
    [[maybe_unused]] double median(t_treeNode<T> *node);
    [[maybe_unused]] double average(t_treeNode<T> *node);
    [[maybe_unused]] double standardDeviation(t_treeNode<T> *node);
    [[maybe_unused]] double standardError(t_treeNode<T> *node);
    [[maybe_unused]] double variance(t_treeNode<T> *node, const string& type);
    [[maybe_unused]] double covariance(t_treeNode<T> *root_, t_treeNode<T> *node, const string& type);
    [[maybe_unused]] double sumOfSquaredMeanDifferences(t_treeNode<T> *node);
    [[maybe_unused]] double skewness(t_treeNode<T> *node);
    [[maybe_unused]] double kurtosis(t_treeNode<T> *node);
    [[maybe_unused]] T mode([[maybe_unused]] t_treeNode<T> *node);
    [[maybe_unused]] double correlationCoefficient(t_treeNode<T> *root_, t_treeNode<T> *node);
    [[maybe_unused]] double zMultiplier(double alpha);
    [[maybe_unused]] pair<double, double> confidenceInterval(t_treeNode<T> *node, double alpha);
    // private destroy methods
    [[maybe_unused]] void destroyTree(t_treeNode<T> *node);
public:
    // constructors
    t_BST();
    t_BST(const t_BST<T> &other);
    t_BST(t_BST<T> &&other) noexcept ;
    t_BST<T> &operator=(const t_BST<T> &other);
    t_BST<T> &operator=(t_BST<T> &&other) noexcept;
    [[maybe_unused]] explicit t_BST([[maybe_unused]] vector<T> &data);
    [[maybe_unused]] t_BST(const string& fileName, char delimiter);

    [[maybe_unused]] t_BST(int totalElements, T min, T max);
    // public insert methods
    [[maybe_unused]] void insert(T data);
    // public remove methods
    [[maybe_unused]] void remove(T data);
    [[maybe_unused]] void removeAll(T data);
    // public search methods
    [[maybe_unused]] t_treeNode<T> *find(T data);
    [[maybe_unused]] t_treeNode<T> *findMax();
    [[maybe_unused]] t_treeNode<T> *findMin();
    [[maybe_unused]] t_treeNode<T> *findParent(T data);
    [[maybe_unused]] t_treeNode<T> *findSuccessor(T data);
    [[maybe_unused]] t_treeNode<T> *findPredecessor(T data);
    T max();
    T min();
    // public print methods
    [[maybe_unused]] void printTree();
    [[maybe_unused]] void inOrder();
    [[maybe_unused]] void preOrder();
    [[maybe_unused]] void postOrder();
    [[maybe_unused]] void levelOrder();
    [[maybe_unused]] void print2D();
    // public getter methods
    [[maybe_unused]] t_treeNode<T> *getRoot();
    [[maybe_unused]] int getSize();
    [[maybe_unused]] int getHeight();
    [[maybe_unused]] int getDiameter();
    [[maybe_unused]] int getBalanceFactor();
    [[maybe_unused]] int getCount(T key);
    [[maybe_unused]] int getCount();
    [[maybe_unused]] int countLeaves();
    // public balance method
    [[maybe_unused]] void balance();
    // public transformation methods
    [[maybe_unused]] void replaceWithSums();
    // public statistics methods
    [[maybe_unused]] double sum();
    [[maybe_unused]] double median();
    [[maybe_unused]] double mean();
    [[maybe_unused]] double average();
    [[maybe_unused]] double standardDeviation();
    [[maybe_unused]] double variance(const string& type = "population");
    [[maybe_unused]] double covariance(t_treeNode<T> *node, string type = "population");
    [[maybe_unused]] double standardError();

    [[maybe_unused]] [[maybe_unused]] double sumOfSquaredMeanDifferences();
    [[maybe_unused]] T mode();
    [[maybe_unused]] double skewness();
    [[maybe_unused]]  double kurtosis();
    [[maybe_unused]] double range();
    [[maybe_unused]] double correlationCoefficient(t_treeNode<T> *node);
    [[maybe_unused]] pair<double, double> confidenceInterval(double confidence);
    // public boolean methods
    [[maybe_unused]] bool isEmpty();
    // public tree filler methods
    [[maybe_unused]] void fillTreeFromFile([[maybe_unused]] const string& filename,
                                           [[maybe_unused]] char delimiter);
    // public destructor methods
    [[maybe_unused]] void destroy();
    ~t_BST();

};
#endif //STATBST_T_BST_H
template <typename T>
[[maybe_unused]] void t_BST<T>::insert(t_treeNode<T> *&root_, t_treeNode<T> *&newNode) {
    if (root_ == nullptr) {
        root_ = newNode;
        size++;
        return;
    }
    if (newNode->data < root_->data) {
        insert(root_->left, newNode);
    } else {
        insert(root_->right, newNode);
    }
    balance(root_);
}
// helper function that will remove a node from the tree with the given data
template <typename T>
[[maybe_unused]] void t_BST<T>::remove(t_treeNode<T> *&root_, T data) {
    if (root_ == nullptr) {
        return;
    }
    if (data < root_->data) {
        remove(root_->left, data);
    } else if (data > root_->data) {
        remove(root_->right, data);
    } else {
        if (root_->left == nullptr && root_->right == nullptr) {
            delete root_;
            root_ = nullptr;
            size--;
        } else if (root_->left == nullptr) {
            t_treeNode<T> *temp = root_;
            root_ = root_->right;
            delete temp;
            size--;
        } else if (root_->right == nullptr) {
            t_treeNode<T> *temp = root_;
            root_ = root_->left;
            delete temp;
            size--;
        } else {
            t_treeNode<T> *temp = findMin(root_->right);
            root_->data = temp->data;
            remove(root_->right, temp->data);
        }
    }
    balance(root_);
}
// helper function that removes all nodes with the given data
template <typename T>
[[maybe_unused]] void t_BST<T>::removeAll(t_treeNode<T> *&root_, T data) {
    if (root_ == nullptr) {
        return;
    }
    do {
        remove(root_, data);
    }while(getCount(data) > 0);
}
// helper function that prints the tree in order
template <typename T>
[[maybe_unused]] void t_BST<T>::printTree(t_treeNode<T> *root_) {
    if (root_ != nullptr) {
        printTree(root_->left);
        cout << root_->data << " ";
        printTree(root_->right);
    }
}
// helper function to print the tree in order
template <typename T>
[[maybe_unused]] void t_BST<T>::inOrder(t_treeNode<T> *root_) {
    if (root_ != nullptr) {
        inOrder(root_->left);
        cout << root_->data << " ";
        inOrder(root_->right);
    }
}
// helper function to print in pre-order
template <typename T>
[[maybe_unused]] void t_BST<T>::preOrder(t_treeNode<T> *root_) {
    if (root_ != nullptr) {
        cout << root_->data << " ";
        preOrder(root_->left);
        preOrder(root_->right);
    }
}
template <typename T>
[[maybe_unused]] void t_BST<T>::postOrder(t_treeNode<T> *node) {
    if (node != nullptr) {
        postOrder(node->left);
        postOrder(node->right);
        cout << node->data << " ";
    }
}
template <typename T>
[[maybe_unused]] void t_BST<T>::levelOrder(t_treeNode<T> *node) {
    if (node == nullptr) {
        return;
    }
    queue<t_treeNode<T> *> q;
    q.push(node);
    while (!q.empty()) {
        t_treeNode<T> *temp = q.front();
        q.pop();
        cout << temp->data << " ";
        if (temp->left != nullptr) {
            q.push(temp->left);
        }
        if (temp->right != nullptr) {
            q.push(temp->right);
        }
    }
}

template <typename T>
[[maybe_unused]] void t_BST<T>::print2DUtil(t_treeNode<T> *root_, int space) {
    if (root_ == nullptr) {
        return;
    }
    space += 10;
    print2DUtil(root_->right, space);
    cout << endl;
    for (int i = 10; i < space; i++) {
        cout << " ";
    }
    cout << root_->data << endl;
    print2DUtil(root_->left, space);
}


template <typename T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::findSuccessor(t_treeNode<T> *node) {
    if (node == nullptr) {
        return nullptr;
    }
    if (node->right != nullptr) {
        return findMin(node->right);
    }
    t_treeNode<T> *parent = findParent(root, node->data);
    while (parent != nullptr && parent->right != nullptr &&
           parent->right->data != node->data) {
        node = parent;
        parent = findParent(root, node->data);
    }
    return parent;
}

template <typename T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::findPredecessor(t_treeNode<T> *node) {
    if (node == nullptr) {
        return nullptr;
    }
    if (node->left != nullptr) {
        return findMax(node->left);
    }
    t_treeNode<T> *parent = findParent(root, node->data);
    while (parent != nullptr && parent->left != nullptr &&
           parent->left->data != node->data) {
        node = parent;
        parent = findParent(root, node->data);
    }
    return parent;
}

// private balance methods
// right rotate
template <typename T>
[[maybe_unused]] void t_BST<T>::rightRotate(t_treeNode<T> *&node) {
    t_treeNode<T> *temp = node->left;
    node->left = temp->right;
    temp->right = node;
    node = temp;
}
// left rotate
template <typename T>
[[maybe_unused]] void t_BST<T>::leftRotate(t_treeNode<T> *&node) {
    t_treeNode<T> *temp = node->right;
    node->right = temp->left;
    temp->left = node;
    node = temp;
}
// rotate left right
template <typename T>
[[maybe_unused]] void t_BST<T>::leftRightRotate(t_treeNode<T> *&node) {
    leftRotate(node->left);
    rightRotate(node);
}
// rotate right left
template <typename T>
[[maybe_unused]] void t_BST<T>::rightLeftRotate(t_treeNode<T> *&node) {
    rightRotate(node->right);
    leftRotate(node);
}
// balance tree
template <typename T>
[[maybe_unused]] void t_BST<T>::balance(t_treeNode<T> *&node) {
    if (node == nullptr) {
        return;
    }
    int balanceFactor = getBalanceFactor(node);
    if (balanceFactor > 1) {
        if (getBalanceFactor(node->left) < 0) {
            leftRotate(node->left);
        }
        rightRotate(node);
    } else if (balanceFactor < -1) {
        if (getBalanceFactor(node->right) > 0) {
            rightRotate(node->right);
        }
        leftRotate(node);
    }
}
// get leaves
template <typename T>
[[maybe_unused]] int t_BST<T>::getLeaves(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    if (node->left == nullptr && node->right == nullptr) {
        return 1;
    }
    return getLeaves(node->left) + getLeaves(node->right);
}
template <typename T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::find(t_treeNode<T> *node, T data) {
    if (node == nullptr) {
        return nullptr;
    }
    if (data == node->data) {
        return node;
    } else if (data < node->data) {
        return find(node->left, data);
    } else {
        return find(node->right, data);
    }
}
template <typename T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::findMin(t_treeNode<T> *node) {
    if (node == nullptr) {
        return nullptr;
    }
    while (node->left != nullptr) {
        node = node->left;
    }
    return node;
}
template <typename T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::findMax(t_treeNode<T> *node) {
    if (node == nullptr) {
        return nullptr;
    }
    while (node->right != nullptr) {
        node = node->right;
    }
    return node;
}
template <typename T>
[[maybe_unused]] void t_BST<T>::destroyTree(t_treeNode<T> *node) {
    if (node != nullptr) {
        destroyTree(node->left);
        destroyTree(node->right);
        delete node;
    }
}
template <typename T>
[[maybe_unused]] int t_BST<T>::getBalanceFactor(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    return getHeight(node->left) - getHeight(node->right);
}
template <typename T>
[[maybe_unused]] int t_BST<T>::getHeight(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    int leftHeight = getHeight(node->left);
    int rightHeight = getHeight(node->right);
    return (leftHeight > rightHeight ? leftHeight : rightHeight) + 1;
}
template <typename T>
[[maybe_unused]] int t_BST<T>::getDiameter(t_treeNode<T>* root_) {
    if (root_ == nullptr) {
        return 0;
    }
    int leftHeight = getHeight(root_->left);
    int rightHeight = getHeight(root_->right);
    int leftDiameter = getDiameter(root_->left);
    int rightDiameter = getDiameter(root_->right);
    return std::max(leftHeight + rightHeight + 1, std::max(leftDiameter, rightDiameter));
}
template <typename T>
[[maybe_unused]] int t_BST<T>::getCount(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    return getCount(node->left) + getCount(node->right) + 1;
}
template <typename T>
[[maybe_unused]] T t_BST<T>::replaceWithSums(t_treeNode<T> *&node) {
    if (node == nullptr) {
        return 0;
    }
    if (node->left == nullptr && node->right == nullptr) {
        return node->data;
    }

    T leftSum = replaceWithSums(node->left);
    T rightSum = replaceWithSums(node->right);
    T temp = node->data;
    node->data = leftSum + rightSum;
    return node->data + temp;
}
template <typename T>
[[maybe_unused]] double t_BST<T>::sum(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    if (node->left == nullptr && node->right == nullptr) {
        return node->data;
    }

    double leftSum = sum(node->left);
    double rightSum = sum(node->right);
    return node->data + leftSum + rightSum;
}
template <typename T>
[[maybe_unused]] double t_BST<T>::median(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    int count = getCount(node);
    int currCount = 0;
    auto *current = node;
    double median = 0;
    t_treeNode<T> *pre;
    while (current != nullptr) {
        if (current->left == nullptr) {
            if (currCount == count / 2) {
                median = current->data;
            }
            currCount++;
            current = current->right;
        } else {
            pre = current->left;
            while (pre->right != nullptr && pre->right != current) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = current;
                current = current->left;
            } else {
                pre->right = nullptr;
                if (currCount == count / 2) {
                    median = current->data;
                }
                currCount++;
                current = current->right;
            }
        }
    }
    return median;
}
template <typename T>
[[maybe_unused]] double t_BST<T>::average(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    int count = getCount(node);
    double sum = this->sum(node);
    return sum / count;
}
template <typename T>
[[maybe_unused]] double t_BST<T>::standardDeviation(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    int count = getCount(node);
    double average = this->average();
    double sum = 0;
    auto *current = node;
    while (current != nullptr) {
        if (current->left == nullptr) {
            sum += pow(current->data - average, 2);
            current = current->right;
        } else {
            auto *pre = current->left;
            while (pre->right != nullptr && pre->right != current) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = current;
                current = current->left;
            } else {
                pre->right = nullptr;
                sum += pow(current->data - average, 2);
                current = current->right;
            }
        }
    }
    return sqrt(sum / count);
}
template <typename T>
[[maybe_unused]] double t_BST<T>::standardError(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    int count = getCount(node);
    double standardDeviation = this->standardDeviation();
    return standardDeviation / sqrt(count);
}
template <typename T>
[[maybe_unused]] double t_BST<T>::variance(t_treeNode<T> *node, const string& type) {
    if (node == nullptr) {
        return 0;
    }
    int count = getCount(node);
    double average = this->average();
    double sum = 0;
    auto *current = node;
    while (current != nullptr) {
        if (current->left == nullptr) {
            sum += pow(current->data - average, 2);
            current = current->right;
        } else {
            auto *pre = current->left;
            while (pre->right != nullptr && pre->right != current) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = current;
                current = current->left;
            } else {
                pre->right = nullptr;
                sum += pow(current->data - average, 2);
                current = current->right;
            }
        }
    }
    if (type == "population") {
        return sum / count;
    } else {
        return sum / (count - 1);
    }
}
template <typename T>
[[maybe_unused]] double t_BST<T>::covariance(t_treeNode<T> *root_, t_treeNode<T> *node, const string& type) {
    if (root_ == nullptr || node == nullptr) {
        return 0;
    }
    int count = getCount(root_);
    double sum_XY = 0;
    double rootAvg = average(root);
    double nodeAvg = average(node);
    vector<T> x;
    vector<T> y;
    // adding the root data to the vector x
    while (root_ != nullptr) {
        if (root_->left == nullptr) {
            x.push_back(root_->data);
            root_ = root_->right;
        } else {
            auto *pre = root_->left;
            while (pre->right != nullptr && pre->right != root_) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = root_;
                root_ = root_->left;
            } else {
                pre->right = nullptr;
                x.push_back(root_->data);
                root_ = root_->right;
            }
        }
    }
    // adding the node data to the vector y
    while (node != nullptr) {
        if (node->left == nullptr) {
            y.push_back(node->data);
            node = node->right;
        } else {
            auto *pre = node->left;
            while (pre->right != nullptr && pre->right != node) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = node;
                node = node->left;
            } else {
                pre->right = nullptr;
                y.push_back(node->data);
                node = node->right;
            }
        }
    }
    // sort the vectors
    sort(x.begin(), x.end());
    sort(y.begin(), y.end());
    // calculate the sum of the x-xAvg * y-yAvg
    for (int i = 0; i < x.size(); i++) {
        sum_XY += (x[i] - rootAvg) * (y[i] - nodeAvg);
    }
    if (type == "population") {
        return sum_XY / count;
    } else {
        return sum_XY / (count - 1);
    }
}

template <typename T>
[[maybe_unused]] double t_BST<T>::sumOfSquaredMeanDifferences(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    double average = this->average(node);
    double sum = 0;
    auto *current = node;
    while (current != nullptr) {
        if (current->left == nullptr) {
            sum += pow(current->data - average, 2);
            current = current->right;
        } else {
            auto *pre = current->left;
            while (pre->right != nullptr && pre->right != current) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = current;
                current = current->left;
            } else {
                pre->right = nullptr;
                sum += pow(current->data - average, 2);
                current = current->right;
            }
        }
    }
    return sum;
}

template <typename T>
[[maybe_unused]] double t_BST<T>::skewness(t_treeNode<T> *node) {
    if (node == nullptr) {
        return 0;
    }
    int count = getCount(node);
    double average = this->average();
    double sum = 0;
    auto *current = node;
    while (current != nullptr) {
        if (current->left == nullptr) {
            sum += pow(current->data - average, 3);
            current = current->right;
        } else {
            auto *pre = current->left;
            while (pre->right != nullptr && pre->right != current) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = current;
                current = current->left;
            } else {
                pre->right = nullptr;
                sum += pow(current->data - average, 3);
                current = current->right;
            }
        }
    }
    return sum / (count * pow(this->standardDeviation(), 3));
}
template <typename T>
[[maybe_unused]] double t_BST<T>::kurtosis(t_treeNode<T> *node) {
    // the central moment of the distribution
    if (node == nullptr) {
        return 0;
    }
    int count = getCount(node);
    double average = this->average();
    double sum = 0;
    auto *current = node;
    while (current != nullptr) {
        if (current->left == nullptr) {
            sum += pow(current->data - average, 4);
            current = current->right;
        } else {
            auto *pre = current->left;
            while (pre->right != nullptr && pre->right != current) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = current;
                current = current->left;
            } else {
                pre->right = nullptr;
                sum += pow(current->data - average, 4);
                current = current->right;
            }
        }
    }
    return sum / (count * pow(this->standardDeviation(), 4));
}
template <typename T>
[[maybe_unused]] T t_BST<T>::mode([[maybe_unused]] t_treeNode<T> *node) {
    // find the value that occurs most often and return it
    if (node == nullptr) {
        return 0;
    }
    unordered_map<T, int> count;
    queue<t_treeNode<T> *> q;
    q.push(node);
    [[maybe_unused]] int maxFreq = 0;
    while (!q.empty()) {
        [[maybe_unused]] auto *current = q.front();
        q.pop();
        if (current->left == nullptr && current->right == nullptr) {
            count[current->data]++;
            if (count[current->data] > maxFreq) {
                maxFreq = count[current->data];
            }
        } else {
            if (current->left != nullptr) {
                q.push(current->left);
            }
            if (current->right != nullptr) {
                q.push(current->right);
            }
        }
    }
    T mode;
    for ([[maybe_unused]] auto &i : count) {
        if (i.second == maxFreq) {
            mode = i.first;
        }
    }
    return mode;
}
template <typename T>
[[maybe_unused]] double t_BST<T>::correlationCoefficient(t_treeNode<T> *root_, t_treeNode<T> *node) {
    if (root_ == nullptr || node == nullptr) {
        return 0;
    }
    double sum_XY = 0;
    double rootAvg = average(root);
    double nodeAvg = average(node);
    double squareSum_X = this->sumOfSquaredMeanDifferences(root);
    double squareSum_Y = this->sumOfSquaredMeanDifferences(node);
    vector<T> x;
    vector<T> y;
    // adding the root data to the vector x
    while (root_ != nullptr) {
        if (root_->left == nullptr) {
            x.push_back(root_->data);
            root_ = root_->right;
        } else {
            auto *pre = root_->left;
            while (pre->right != nullptr && pre->right != root_) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = root_;
                root_ = root_->left;
            } else {
                pre->right = nullptr;
                x.push_back(root_->data);
                root_ = root_->right;
            }
        }
    }
    // adding the node data to the vector y
    while (node != nullptr) {
        if (node->left == nullptr) {
            y.push_back(node->data);
            node = node->right;
        } else {
            auto *pre = node->left;
            while (pre->right != nullptr && pre->right != node) {
                pre = pre->right;
            }
            if (pre->right == nullptr) {
                pre->right = node;
                node = node->left;
            } else {
                pre->right = nullptr;
                y.push_back(node->data);
                node = node->right;
            }
        }
    }
    // sort the vectors
    sort(x.begin(), x.end());
    sort(y.begin(), y.end());
    // calculate the sum of the x-xAvg * y-yAvg
    for (int i = 0; i < x.size(); i++) {
        sum_XY += (x[i] - rootAvg) * (y[i] - nodeAvg);
    }
    return sum_XY / sqrt(squareSum_X * squareSum_Y);
}
template <typename T>
[[maybe_unused]] double t_BST<T>::zMultiplier(double alpha) {
    // calculate the z multiplier for the given alpha
    map<double, double> z;
    //read the values from the file where each line is a pair of values
    // separated by a space
    ifstream file("zscores.txt");
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        double x, y;
        ss >> x >> y;
        z[x] = y;
    }
    return z[alpha];
}

template <typename T>
[[maybe_unused]] pair<double, double> t_BST<T>::confidenceInterval(t_treeNode<T> *node, double alpha) {
    // return the confidence interval of the data in the node
    // alpha is the significance level
    double average = this->average(node);
    double standardDeviation = this->standardDeviation(node);
    double z = this->zMultiplier(alpha);
    double confidenceInterval = z * standardDeviation / sqrt(getCount(node));
    double lowerBound = average - confidenceInterval;
    double upperBound = average + confidenceInterval;
    return make_pair(lowerBound, upperBound);
}

template <typename T>
[[maybe_unused]] int t_BST<T>::keyCountUtil(t_treeNode<T>* root_, T key) {
    if (root_ == nullptr) {
        return 0;
    }
    int count = 0;
    if (root_->data == key) {
        count++;
    }
    count += keyCountUtil(root_->left, key);
    count += keyCountUtil(root_->right, key);
    return count;
}

//____________________________public functions__________________________________
// constructor
template <typename T>
[[maybe_unused]] t_BST<T>::t_BST() {
    root = nullptr;
    size = 0;
}
// copy constructor
template<class T>
[[maybe_unused]] t_BST<T>::t_BST(const t_BST<T> &other) {
    root = copy(other.root);
}
// move constructor
template<class T>
[[maybe_unused]] t_BST<T>::t_BST(t_BST<T> &&other) noexcept {
    root = other.root;
    other.root = nullptr;
}
// copy assignment operator
template<class T>
[[maybe_unused]] t_BST<T> &t_BST<T>::operator=(const t_BST<T> &other) {
    if (this != &other) {
        destroy();
        copyTree(other.root, root);
    }
    return *this;
}
// move assignment operator
template<class T>
[[maybe_unused]] t_BST<T> &t_BST<T>::operator=(t_BST<T> &&other) noexcept {
    if (this != &other) {
        destroy();
        root = other.root;
        other.root = nullptr;
    }
    return *this;
}
// constructor with initializer list
template<class T>
[[maybe_unused]] t_BST<T>::t_BST([[maybe_unused]] vector<T> &data) {
    root = nullptr;
    size = 0;
    for ([[maybe_unused]] auto &i : data) {
        insert(i);
    }
}
// constructor to create a BST from data in a file
template<class T>
[[maybe_unused]] t_BST<T>::t_BST(const string& fileName, char delimiter) {
    root = nullptr;
    size = 0;
    this->fillTreeFromFile(fileName, delimiter);
}

template<class T>
[[maybe_unused]] t_BST<T>::t_BST(int totalElements, T min, T max) {
    // fill the tree with random numbers in the range [min, max]
    root = nullptr;
    size = 0;
    // seed the random number generator from the system clock and random library
    srandom(time(nullptr));
    for (int i = 0; i < totalElements; i++) {
        insert(random() % (max - min + 1) + min);
    }
}
/**
 * @brief insert data into the tree
 * @param data  data to be inserted
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::insert(T data) {
    auto *node = new t_treeNode<T>(data);
    insert(root, node);
}
/**
 * @brief destroy tree
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::destroy() {
    destroyTree(root);
    root = nullptr;
    size = 0;
}
/**
* @brief find the node with the given data
* @param data  the data to find
* @return  the node with the given data
*/
template <typename T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::find(T data) {
    return find(root, data);
}
/**
 * @brief find the maximum node
 * @return  the maximum node
 */
template <typename T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::findMax() {
    return findMax(root);
}
template <typename T>
[[maybe_unused]] T t_BST<T>::max() {
    return findMax()->data;
}
/**
 * @brief find the minimum node
 * @return  the minimum node
 */
template <typename T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::findMin() {
    return findMin(root);
}
template <typename T>
[[maybe_unused]] T t_BST<T>::min() {
    return findMin()->data;
}
/**
* @brief finds the parent of the node, if the node is not in the tree,
* returns nullptr
* @param data  the data of the node to find the parent of
* @return  the parent of the node, if the node is not in the tree,
*  returns nullptr
*/
template<class T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::findParent(T data) {
    t_treeNode<T> *node = find(data);
    if (node == nullptr) {
        return nullptr;
    }
    if (node->parent == nullptr) {
        return nullptr;
    }
    return node->parent;
}
/**
 * @brief finds the successor of a node, or the leftmost node of the right subtree
 * @param data  the data of the node to find the successor of
 * @return the successor of the node
 */
template<class T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::findSuccessor(T data) {
    t_treeNode<T> *node = find(data);
    if (node == nullptr) {
        return nullptr;
    }
    if (node->right != nullptr) {
        return findMin(node->right);
    }
    t_treeNode<T> *parent = findParent(data);
    while (parent != nullptr && parent->data < data) {
        node = parent;
        parent = findParent(data);
    }
    return node;
}
/**
 * @brief finds the predecessor of a node, or the right most child of the
 * left subtree of the node
 * @param data  the data of the node to find the predecessor of
 * @return the predecessor of the node
 */
template<class T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::findPredecessor(T data) {
    t_treeNode<T> *node = find(data);
    if (node == nullptr) {
        return nullptr;
    }
    if (node->left != nullptr) {
        return findMax(node->left);
    }
    t_treeNode<T> *parent = findParent(data);
    while (parent != nullptr && parent->left != node) {
        node = parent;
        parent = findParent(data);
    }
    return parent;
}
/**
 * @brief counts how many times a key appears in the tree
 * @param key  the key to count
 * @return  the number of times the key appears in the tree
 */
template <typename T>
[[maybe_unused]] int t_BST<T>::getCount(T key) {
    return keyCountUtil(root, key);
}

/**
 * @brief remove a node from the tree
 * @param data  the data of the node to remove
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::remove(T data) {
    remove(root, data);
}
template <typename T>
[[maybe_unused]] void t_BST<T>::removeAll(T data) {
    removeAll(root, data);
}

/**
 *  @brief prints the tree in order
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::printTree() {
    printTree(root);
}
/**
 *  @brief print the tree in order
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::inOrder() {
    inOrder(root);
}
/**
 *  @brief print the tree in pre order
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::preOrder() {
    preOrder(root);
}
/**
 *  @brief print the tree in post order
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::postOrder() {
    postOrder(root);
}
/**
 *  @brief print the tree in level order
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::levelOrder() {
    levelOrder(root);
}
/**
 *  @brief prints the tree in a 2d fashion
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::print2D() {
    print2DUtil(root, 0);
}

// getters
/**
 * @brief get the root of the tree
 *  @return the root of the tree
 */
template <typename T>
[[maybe_unused]] t_treeNode<T> *t_BST<T>::getRoot() {
    return root;
}
/**
 * @brief get the size of the tree
 * @return the size of the tree
 */
template <typename T>
[[maybe_unused]] int t_BST<T>::getSize() {
    return size;
}
/**
 *  @brief get the height of the tree
 * @return  the height of the tree from the node
 */
template <typename T>
[[maybe_unused]] int t_BST<T>::getHeight() {
    return getHeight(root);
}

/**
 * @brief get the diameter of the tree
 * @return  the diameter of the tree
 */
template <typename T>
[[maybe_unused]] int t_BST<T>::getDiameter() {
    return getDiameter(root);
}

/**
 *  @brief get the balance factor of the tree
 * @return  the balance factor of the tree from the node
 */
template <typename T>
[[maybe_unused]] int t_BST<T>::getBalanceFactor() {
    return getBalanceFactor(root);
}

/**
 *  @brief balance the tree
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::balance() {
    balance(root);
}
/**
 * @brief replaces each node with the sum of its left and right subtrees
 */
template <typename T>
[[maybe_unused]] void t_BST<T>::replaceWithSums() {
    replaceWithSums(root);
}
/**
 * @brief counts the number of leaves in the tree
 * @return  the number of leaves in the tree
 */
template <typename T>
[[maybe_unused]] int t_BST<T>::countLeaves() {
    return getLeaves(root);
}
template <typename T>
[[maybe_unused]] int t_BST<T>::getCount() {
    return getCount(root);
}

/**
 * @brief sums the nodes in the tree
 * @return  the sum of the nodes in the tree
 */
template <typename T>
[[maybe_unused]] double t_BST<T>::sum() {
    return sum(root);
}

// median
/**
 * @brief finds the median of the tree
 * @return  the median of the tree
 */
template <typename T>
[[maybe_unused]] double t_BST<T>::median() {
    return median(root);
}
template <typename T>
[[maybe_unused]] double t_BST<T>::mean() {
    return this->average();
}
/**
 * @brief finds the average of the tree
 * @return  the average of the tree
 */
template <typename T>
[[maybe_unused]] double t_BST<T>::average() {
    return average(root);
}
/**
 * @brief finds the standard deviation of the tree
 * @return  the standard deviation of the tree
 */
template <typename T>
[[maybe_unused]] double t_BST<T>::standardDeviation() {
    return standardDeviation(root);
}
template <typename T>
[[maybe_unused]] double t_BST<T>::variance(const string& type) {
    return variance(root, type);
}
template <typename T>
[[maybe_unused]] double t_BST<T>::covariance(t_treeNode<T> *node, string type) {
    return covariance(root, node, type);
}

/**
 * @brief finds the standard error of the tree
 * @return  the standard error of the tree
 */
template <typename T>
[[maybe_unused]] double t_BST<T>::standardError() {
    return standardError(root);
}
template <typename T>
[[maybe_unused]] double t_BST<T>::sumOfSquaredMeanDifferences() {
    return sumOfSquaredMeanDifferences(root);
}

/**
 * @brief finds the mode of the tree
 * @return  the mode of the tree
 */
template <typename T>
[[maybe_unused]] T t_BST<T>::mode() {
    return mode(root);
}
/**
 * @brief finds the skewness, or the ratio of the difference between the
 * mean and the median to the standard deviation
 * @return
 */
template <typename T>
[[maybe_unused]] double t_BST<T>::skewness() {
    return skewness(root);
}
template <typename T>
[[maybe_unused]] double t_BST<T>::kurtosis() {
    return kurtosis(root);
}
template <typename T>
[[maybe_unused]] double t_BST<T>::range() {
    // difference between max and min
    return max() - min();
}
/**
 * @brief finds the correlation coefficient between two trees
 * @param node  the node to compare to
 * @return  the correlation coefficient between the two trees
 */
template <typename T>
[[maybe_unused]] double t_BST<T>::correlationCoefficient(t_treeNode<T> *node) {
    return correlationCoefficient(root, node);
}

/**
 * @brief finds the confidence interval of the tree
 * @param confidence  the confidence level as a decimal
 * @return  a pair of doubles that represent the lower and upper bounds
 * of the confidence interval
 */
template <typename T>
[[maybe_unused]] pair<double, double> t_BST<T>::confidenceInterval(double confidence) {
    return confidenceInterval(root, confidence);
}
/**
 * @brief checks to see if the tree is empty
 * @return true if tree is empty else false
 */
template <typename T>
[[maybe_unused]] bool t_BST<T>::isEmpty() {
    return root == nullptr;
}
template <typename T>
[[maybe_unused]] void t_BST<T>::fillTreeFromFile([[maybe_unused]] const string& filename,
                                                 [[maybe_unused]] char delimiter) {
    ifstream file;
    file.open(filename);
    string line;
    vector<T> tokens;
    while (getline(file, line)) {
        stringstream ss(line);
        string token;
        while (getline(ss, token, delimiter)) {
            tokens.push_back(stod(token));
        }
    }
    file.close();
    for ([[maybe_unused]] auto& token : tokens) {
        insert(token);
    }
}
/**
 * destructor
 */
template <typename T>
t_BST<T>::~t_BST() {
    destroyTree(root);
}

