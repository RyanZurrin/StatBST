#include "T_BST.h"

int main() {
    // test the BST
    t_BST<int> bst(30, 0, 50);
    // fill bst from file named infile.txt with a comma delimited file
    // containing integers
    //bst.fillTreeFromFile("infile.txt", ',');
    t_BST<int> bst2(30, 0, 50);

    // fill bst2 from file named infile2.txt with a space delimited file
    // containing integers
    //bst2.fillTreeFromFile("infile2.txt", ' ');
    // getCount how many times 10 is in the BST
    cout << "10 appears " << bst.getCount(10) << " times" << endl;
    cout << "10 appears " << bst2.getCount(10) << " times" << endl;
    // get the size of the BST
    cout << "The BST has " << bst.getSize() << " nodes" << endl;
    cout << "The BST2 has " << bst2.getSize() << " nodes" << endl;
    // get the diameter of the BST
    cout << "The diameter of the BST is " << bst.getDiameter() << endl;
    cout << "The diameter of the BST2 is " << bst2.getDiameter() << endl;
    // get the sum of tree
    cout << "The sum of the BST is " << bst.sum() << endl;
    cout << "The sum of the BST2 is " << bst2.sum() << endl;
    // get the median of the BST
    cout << "The median of the BST is " << bst.median() << endl;
    cout << "The median of the BST2 is " << bst2.median() << endl;
    // average of the BST
    cout << "The average of the BST is " << bst.average() << endl;
    cout << "The average of the BST2 is " << bst2.average() << endl;
    // std deviation of the BST
    cout << "The standard deviation of the BST is " <<
    bst.standardDeviation() << endl;
    cout << "The standard deviation of the BST2 is " <<
    bst2.standardDeviation() << endl;
    // standard error of the BST
    cout << "The standard error of the BST is " << bst.standardError() << endl;
    cout << "The standard error of the BST2 is " << bst2.standardError() << endl;
    // get the mode of the BST
    cout << "The mode of the BST is " << bst.mode() << endl;
    cout << "The mode of the BST2 is " << bst2.mode() << endl;
    // get the range of the BST
    cout << "The range of the BST is " << bst.range() << endl;
    cout << "The range of the BST2 is " << bst2.range() << endl;
    // get the skewness of the BST
    cout << "The skewness of the BST is " << bst.skewness() << endl;
    cout << "The skewness of the BST2 is " << bst2.skewness() << endl;
    // get the kurtosis of the BST
    cout << "The kurtosis of the BST is " << bst.kurtosis() << endl;
    cout << "The kurtosis of the BST2 is " << bst2.kurtosis() << endl;
    // get the sum of squared mean differences
    cout << "The sum of squared mean differences is " <<
    bst.sumOfSquaredMeanDifferences() << endl;
    cout << "The sum of squared mean differences is " <<
    bst2.sumOfSquaredMeanDifferences() << endl;
    // get confidence interval for 95% confidence
    pair<double, double> ci = bst.confidenceInterval(0.95);
    pair<double, double> ci2 = bst2.confidenceInterval(0.95);
    // covariance of the BST
    cout << "The covariance of the BST and BST2 is " <<
    bst.covariance(bst2.getRoot()) << endl;
    // covariance as sample
    cout << "The covariance sample of the BST and BST2 is " <<
    bst.covariance(bst2.getRoot(),"sample") << endl;
    // variance of the BST
    cout << "The variance of the BST is " << bst.variance("sample") << endl;
    cout << "The variance of the BST2 is " << bst2.variance() << endl;

    cout << "The 95% confidence interval is [" <<
    ci.first << ", " << ci.second << "]" << endl;
    cout << "The 95% confidence interval is [" <<
    ci2.first << ", " << ci2.second << "]" << endl;
    // find confidence interval for 99% confidence and 80 using new pairs
    pair<double, double> ci3 = bst.confidenceInterval(0.99);
    pair<double, double> ci4 = bst2.confidenceInterval(0.99);
    pair<double, double> ci5 = bst.confidenceInterval(0.8);
    pair<double, double> ci6 = bst2.confidenceInterval(0.8);
    cout << "The 99% confidence interval is [" <<
    ci3.first << ", " << ci3.second << "]" << endl;
    cout << "The 99% confidence interval is [" <<
    ci4.first << ", " << ci4.second << "]" << endl;
    cout << "The 80% confidence interval is [" <<
    ci5.first << ", " << ci5.second << "]" << endl;
    cout << "The 80% confidence interval is [" <<
    ci6.first << ", " << ci6.second << "]" << endl;
    // get the correlation of the BST
    cout << "The correlation of the BST is " << bst.correlationCoefficient(
            bst2.getRoot()) << endl;
    // find the successor of the BST
    t_treeNode<int>* t = bst.findMin();
    t_treeNode<int>* test = bst.findInorderSuccessor(t);
    cout << "The successor of " << t->data <<" is " << test->data << endl;


    cout << "BST: " << endl;
    bst.printTree();
    cout << endl;
    cout << "BST2: " << endl;
    bst2.printTree();
    //print the balance factor
    cout << "\nBalance factor: " << bst.getBalanceFactor() << endl;

    // print level order traversal
    cout << "Level order traversal: " << endl;
    bst.levelOrder();

    // print 2d traversal
//    cout << "\n2d traversal: " << endl;
//    bst.print2D();
    // print the height of the tree
    cout << "Height of the tree: " << bst.getHeight() << endl;

    // delete the node with value 10 and print the tree
    bst.removeAll(10);
    cout << "After deleting 10: " << endl;
    bst.printTree();
    // find min and max
    cout << "\nMin: " << bst.min() << endl;
    cout << "Max: " << bst.max() << endl;
    // replace with sums of all nodes
    bst.replaceWithSums();
    cout << "After replacing with sum: " << endl;
    // level order traversal
    bst.levelOrder();



    return 0;
}