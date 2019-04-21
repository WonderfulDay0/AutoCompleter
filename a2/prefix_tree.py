"""CSC148 Assignment 2: Autocompleter classes

=== CSC148 Fall 2018 ===
Department of Computer Science,
University of Toronto

=== Module Description ===
This file contains the design of a public interface (Autocompleter) and two
implementation of this interface, SimplePrefixTree and CompressedPrefixTree.
You'll complete both of these subclasses over the course of this assignment.

As usual, be sure not to change any parts of the given *public interface* in the
starter code---and this includes the instance attributes, which we will be
testing directly! You may, however, add new private attributes, methods, and
top-level functions to this file.
"""
from __future__ import annotations
from typing import Any, List, Optional, Tuple


################################################################################
# The Autocompleter ADT
################################################################################
class Autocompleter:
    """An abstract class representing the Autocompleter Abstract Data Type.
    """
    def __len__(self) -> int:
        """Return the number of values stored in this Autocompleter."""
        raise NotImplementedError

    def insert(self, value: Any, weight: float, prefix: List) -> None:
        """Insert the given value into this Autocompleter.

        The value is inserted with the given weight, and is associated with
        the prefix sequence <prefix>.

        If the value has already been inserted into this prefix tree
        (compare values using ==), then the given weight should be *added* to
        the existing weight of this value.

        Preconditions:
            weight > 0
            The given value is either:
                1) not in this Autocompleter
                2) was previously inserted with the SAME prefix sequence
        """
        raise NotImplementedError

    def autocomplete(self, prefix: List,
                     limit: Optional[int] = None) -> List[Tuple[Any, float]]:
        """Return up to <limit> matches for the given prefix.

        The return value is a list of tuples (value, weight), and must be
        ordered in non-increasing weight. (You can decide how to break ties.)

        If limit is None, return *every* match for the given prefix.

        Precondition: limit is None or limit > 0.
        """
        raise NotImplementedError

    def remove(self, prefix: List) -> None:
        """Remove all values that match the given prefix.
        """
        raise NotImplementedError


################################################################################
# SimplePrefixTree (Tasks 1-3)
################################################################################
class SimplePrefixTree(Autocompleter):
    """A simple prefix tree.

    This class follows the implementation described on the assignment handout.
    Note that we've made the attributes public because we will be accessing them
    directly for testing purposes.

    === Attributes ===
    value:
        The value stored at the root of this prefix tree, or [] if this
        prefix tree is empty.
    weight:
        The weight of this prefix tree. If this tree is a leaf, this attribute
        stores the weight of the value stored in the leaf. If this tree is
        not a leaf and non-empty, this attribute stores the *aggregate weight*
        of the leaf weights in this tree.
    subtrees:
        A list of subtrees of this prefix tree.
    wgt_type:
        A string that represents how aggregate is calculated
    leafs:
        The number of leaves this subtree has

    === Representation invariants ===
    - self.weight >= 0

    - self.leafs >= 0

    - (EMPTY TREE):
        If self.weight == 0, then self.value == [] and self.subtrees == [].
        This represents an empty simple prefix tree.
    - (LEAF):
        If self.subtrees == [] and self.weight > 0, this tree is a leaf.
        (self.value is a value that was inserted into this tree.)
    - (NON-EMPTY, NON-LEAF):
        If len(self.subtrees) > 0, then self.value is a list (*common prefix*),
        and self.weight > 0 (*aggregate weight*).

    - ("prefixes grow by 1")
      If len(self.subtrees) > 0, and subtree in self.subtrees, and subtree
      is non-empty and not a leaf, then

          subtree.value == self.value + [x], for some element x

    - self.subtrees does not contain any empty prefix trees.
    - self.subtrees is *sorted* in non-increasing order of their weights.
      (You can break ties any way you like.)
      Note that this applies to both leaves and non-leaf subtrees:
      both can appear in the same self.subtrees list, and both have a `weight`
      attribute.
    """
    value: Any
    weight: float
    wgt_type: str
    leafs: int
    subtrees: List[SimplePrefixTree]

    def __init__(self, weight_type: str) -> None:
        """Initialize an empty simple prefix tree.

        Precondition: weight_type == 'sum' or weight_type == 'average'.

        The given <weight_type> value specifies how the aggregate weight
        of non-leaf trees should be calculated (see the assignment handout
        for details).
        """
        self.value = []
        self.weight = 0
        self.subtrees = []
        self.leafs = 0
        self.wgt_type = weight_type

    def __len__(self) -> int:
        """Return the number of values stored in this SimplePrefixTree."""
        if self.is_leaf():
            return 1
        else:
            leaf_count = 0
            for subtree in self.subtrees:
                leaf_count += len(subtree)
            return leaf_count

    def insert(self, value: Any, weight: float, prefix: List) -> None:
        """Insert the given value into this SimplePrefixTree.

        The value is inserted with the given weight, and is associated with
        the prefix sequence <prefix>.

        If the value has already been inserted into this prefix tree
        (compare values using ==), then the given weight should be added to
        the existing weight of this value.

        Preconditions:
            weight > 0
            The given value is either:
                1) not in this SimplePrefixTree
                2) was previously inserted with the SAME prefix sequence
        """

        # Base Case:
        # Prefix sequence being empty (We need to insert)
        if not prefix:

            # Set a flag that detects if the value we are adding is already
            # inside the tree
            found = False
            for subtree in self.subtrees:
                if subtree.value == value:

                    # Value presents in Tree
                    found = True

                    # Update the subtree's weight
                    subtree.weight += weight

            # Value does not present in Tree -> Need to add new value
            # Note this must be added as a leaf.
            if not found:
                # Create the new SimplePrefixTree
                leaf = SimplePrefixTree(self.wgt_type)
                leaf.value = value
                leaf.weight = weight
                leaf.leafs = 1

                # Record the current self.leafs for later calculation

                # Append the leaf to self.subtrees
                self.subtrees.append(leaf)

            # Weight Update
            self.leafs = sum([s.leafs for s in self.subtrees])
            if self.wgt_type == 'sum':
                self.weight = sum([s.weight for s in self.subtrees])
            else:
                self.weight = sum([s.weight * s.leafs for s in self.subtrees]) \
                              / self.leafs

        # Recursion Case (Prefix sequence is not empty)
        else:
            # Create a partial sequence, which consists of self.value and the
            # first index of sequence
            par_seq = self.value + prefix[:1]

            # Create a flag that checks if such partial sequence already
            # presents in the current self.subtrees
            in_tree = False
            for subtree in self.subtrees:
                if subtree.value == par_seq:

                    in_tree = True

                    # Record the subtree.leafs for later calculation
                    subtree.insert(value, weight, prefix[1:])

                    # Update attributes of self
                    self.leafs = sum([sub.leafs for sub in self.subtrees])

                    if self.wgt_type == 'sum':
                        self.weight = sum([s.weight for s in self.subtrees])
                    else:
                        self.weight = sum([sub.weight * sub.leafs for sub in
                                           self.subtrees]) / self.leafs
            # If the partial sequence does not present in tree, we need to add
            # an internal tree.
            if not in_tree:
                internal = SimplePrefixTree(self.wgt_type)
                internal.value = par_seq

                old_leafs = self.leafs

                self.subtrees.append(internal)

                # Insert the value with the rest of the sequence to internal
                internal.insert(value, weight, prefix[1:])

                # Update attributes of self
                self.leafs += 1
                if self.wgt_type == 'sum':
                    self.weight += weight
                else:
                    # Similar calculation as above, note self.leafs and
                    # subtree.leafs are already updated above
                    old_weight = old_leafs * self.weight
                    self.weight = (old_weight + weight) / self.leafs

        self.subtrees = sorted(self.subtrees,
                               key=lambda sub: sub.weight, reverse=True)

    def __contains__(self, item: Any) -> bool:
        if self.is_empty():
            return False
        elif self.value == item:
            return True
        else:
            for subtree in self.subtrees:
                if subtree.__contains__(item):
                    return True
            return False

    def autocomplete(self, prefix: List,
                     limit: Optional[int] = None) -> List[Tuple[Any, float]]:
        """Return up to <limit> matches for the given prefix.

        The return value is a list of tuples (value, weight), and must be
        ordered in non-increasing weight. (You can decide how to break ties.)

        If limit is None, return *every* match for the given prefix.

        Precondition: limit is None or limit > 0.
        """
        values = []
        if prefix == []:
            for subtree in self.subtrees:
                values.extend(subtree._find_leaves())
                if limit is not None and limit <= len(values):
                    break
        else:
            if self.value + prefix[0:1] in self:
                for subtree in self.subtrees:
                    if subtree.value == self.value + [prefix[0]]:
                        values.extend(subtree.autocomplete(prefix[1:], limit))
            else:
                return []

        values = sorted(values,
                        key=lambda sort_wgt: sort_wgt[1], reverse=True)

        if limit is not None:
            return values[0:limit]
        else:
            return values

    def _find_leaves(self) -> List[Tuple[Any, float]]:
        """Returns a two-tuple of all the leafs and their weights in this
        prefix tree.
        """
        values = []
        if self.is_leaf():
            values.append((self.value, self.weight))
        else:
            for subtree in self.subtrees:
                values.extend(subtree._find_leaves())

        return values

    def remove(self, prefix: List) -> None:
        """Remove all values that match the given prefix.
        """

        # Base Case, prefix sequence is empty (Self.leafs handled in recursion)
        if prefix == []:
            self.subtrees = []
            self.weight = 0
            self.leafs = 0

        # Prefix sequence non empty -> Recursion
        else:
            par_seq = self.value + prefix[:1]

            for subtree in self.subtrees:
                if subtree.value == par_seq:
                    subtree.remove(prefix[1:])

                    # Weight Update
                    self.leafs = sum([sub.leafs for sub in self.subtrees])

                    if self.wgt_type == 'sum':
                        self.weight = sum([s.weight for s in self.subtrees])
                    else:
                        self.weight = sum([sub.weight*sub.leafs for sub in
                                           self.subtrees]) / self.leafs

                # Remove empty tree (subtrees == [])
                if not subtree.subtrees and subtree.weight == 0:
                    self.subtrees.remove(subtree)

    def is_empty(self) -> bool:
        """Return whether this simple prefix tree is empty."""
        return self.weight == 0.0

    def is_leaf(self) -> bool:
        """Return whether this simple prefix tree is a leaf."""
        return self.weight > 0 and self.subtrees == []

    def __str__(self) -> str:
        """Return a string representation of this tree.

        You may find this method helpful for debugging.
        """
        return self._str_indented()

    def _str_indented(self, depth: int = 0) -> str:
        """Return an indented string representation of this tree.

        The indentation level is specified by the <depth> parameter.
        """
        if self.is_empty():
            return ''
        else:
            s = '  ' * depth + f'{self.value} ({self.weight})\n'
            for subtree in self.subtrees:
                s += subtree._str_indented(depth + 1)
            return s


################################################################################
# CompressedPrefixTree (Task 6)
################################################################################
class CompressedPrefixTree(SimplePrefixTree):
    """A compressed prefix tree implementation.

    While this class has the same public interface as SimplePrefixTree,
    (including the initializer!) this version follows the implementation
    described on Task 6 of the assignment handout, which reduces the number of
    tree objects used to store values in the tree.

    === Attributes ===
    value:
        The value stored at the root of this prefix tree, or [] if this
        prefix tree is empty.
    weight:
        The weight of this prefix tree. If this tree is a leaf, this attribute
        stores the weight of the value stored in the leaf. If this tree is
        not a leaf and non-empty, this attribute stores the *aggregate weight*
        of the leaf weights in this tree.
    subtrees:
        A list of subtrees of this prefix tree.
    weight_type:
        The weight type of this prefix tree.
    leafs:
        The number of leafs in this prefix tree.

    === Representation invariants ===
    - self.weight >= 0

    - (EMPTY TREE):
        If self.weight == 0, then self.value == [] and self.subtrees == [].
        This represents an empty simple prefix tree.
    - (LEAF):
        If self.subtrees == [] and self.weight > 0, this tree is a leaf.
        (self.value is a value that was inserted into this tree.)
    - (NON-EMPTY, NON-LEAF):
        If len(self.subtrees) > 0, then self.value is a list (*common prefix*),
        and self.weight > 0 (*aggregate weight*).

    - **NEW**
      This tree does not contain any compressible internal values.
      (See the assignment handout for a definition of "compressible".)

    - self.subtrees does not contain any empty prefix trees.
    - self.subtrees is *sorted* in non-increasing order of their weights.
      (You can break ties any way you like.)
      Note that this applies to both leaves and non-leaf subtrees:
      both can appear in the same self.subtrees list, and both have a `weight`
      attribute.
    """
    value: Optional[Any]
    weight: float
    subtrees: List[CompressedPrefixTree]
    weight_type: str
    leafs: int

    def __init__(self, weight_type: str) -> None:
        """Initialize an empty simple prefix tree.

        Precondition: weight_type == 'sum' or weight_type == 'average'.

        The given <weight_type> value specifies how the aggregate weight
        of non-leaf trees should be calculated (see the assignment handout
        for details).


        """
        super().__init__(weight_type)
        self.weight_type = weight_type

    def _pre_in_sub(self, prefix: List) -> bool:
        """Returns True if <self.value> contains the longest <prefix>.
        <prefix> is the longest common prefix.
        """
        for i in range(len(prefix)):
            if i >= len(self.value) or prefix[i] != self.value[i]:
                return False
        return True

    def insert(self, value: Any, weight: float, prefix: List) -> None:
        """Insert the given value into this SimplePrefixTree.

        The value is inserted with the given weight, and is associated with
        the prefix sequence <prefix>.

        If the value has already been inserted into this prefix tree
        (compare values using ==), then the given weight should be added to
        the existing weight of this value.

        Preconditions:
            weight > 0
            The given value is either:
                1) not in this SimplePrefixTree
                2) was previously inserted with the SAME prefix sequence

        """
        # If the prefix is an empty list (prefix == [])
        # If the tree is empty (should only occur once)
        if self.is_empty():
            # Set self to prefix
            self.value = prefix

            # Insert the leaf containing the value into self
            # This should update weight
            # self.insert(value, weight, [])
            self._leaf_gen(value, weight)
        # If self.value is an empty list (self.value == [])
        elif not self.value:
            # <long_com> is the largest common prefix of self.subtrees and
            # prefix
            long_com = self._find_common(prefix)

            # Find which subtree contains long_com (if it's not empty) and
            # recurse into that subtree
            found = False
            for subtree in self.subtrees:
                if long_com and subtree._pre_in_sub(long_com):
                    found = True
                    subtree.insert(value, weight, prefix)
                    break

            # If long_com is an empty list (prefix == [])
            # Need to append (prefix tree --> leaf) directly to self.subtrees
            if not found:
                internal = self._create_int(prefix)

                # Weight of internal should be updated by "if not prefix" case
                # internal.insert(value, weight, [])
                internal._leaf_gen(value, weight)
                self.subtrees.append(internal)

            self.leafs = sum([sub.leafs for sub in self.subtrees])
            # Sum
            self._wgt_upd(weight)
        elif self.value:
            # <long_com_self> is the largest common prefix of self.value and
            # prefix
            long_com_self = self._find_common2(prefix)

            # If both the prefix and the longest common prefix are bigger than
            # self.value, then we know that longest common prefix should be a
            # parent of both value AND self
            if len(self.value) > len(long_com_self) < len(prefix):
                # This is the internal node containing a prefix of both
                internal = self._create_int(long_com_self)
                internal.weight, internal.leafs = self.weight, self.leafs

                # This is the internal node containing prefix of value
                int_of_internal = self._create_int(prefix)

                self.value, internal.value = internal.value, self.value

                for subtree in self.subtrees:
                    internal.subtrees.append(subtree)

                self.subtrees = [internal]
                # Maybe order matters, careful (put this after int_of_internal
                # insert
                self.subtrees.append(int_of_internal)

                # Insert value into <prefix> node
                # This also sets the weight of it
                # int_of_internal.insert(value, weight, [])
                int_of_internal._leaf_gen(value, weight)

                self.leafs += 1
                self._wgt_upd(weight)
            # If the both long_com_self is the prefix, then prefix is the parent
            elif len(self.value) > len(long_com_self) == len(prefix):
                internal = self._create_int(prefix)
                internal.weight, internal.leafs = self.weight, self.leafs

                self.value, internal.value = internal.value, self.value

                for subtree in self.subtrees:
                    internal.subtrees.append(subtree)

                self.subtrees = [internal]

                # Should update self.weight
                # self.insert(value, weight, [])
                self._leaf_gen(value, weight)
            # The length of self.value == length of long_com_self, this implies
            # self.value == prefix, so recurse up to "if not prefix" case
            elif len(self.value) == len(long_com_self) == len(prefix):
                # self.insert(value, weight, [])
                self._leaf_gen(value, weight)
            # The prefix should be inserted AFTER self because its strictly
            # less than long_com_self
            elif len(self.value) < len(prefix) > len(long_com_self):
                self.last_case(prefix, long_com_self, value, weight)

        self.subtrees = sorted(self.subtrees,
                               key=lambda sub: sub.weight, reverse=True)

    def last_case(self, prefix: List, long_com_self: List, value: Any,
                  weight: float) -> None:
        """Helper method to literally make the length of the code shorter in one
        block so that the popcorn PyTA thingy doesn't say it's too long.
        *Represents the last case of insert*.
        """
        found = False
        for subtree in self.subtrees:
            sub_long_com = subtree._find_common2(prefix)
            if len(long_com_self) < len(sub_long_com):
                found = True
                subtree.insert(value, weight, prefix)
                break

        if not found:
            internal = self._create_int(prefix)

            # internal.insert(value, weight, [])
            internal._leaf_gen(value, weight)

            self.subtrees.append(internal)

        self.leafs = sum([sub.leafs for sub in self.subtrees])

        self._wgt_upd(weight)

    def _wgt_upd(self, weight: float) -> None:
        """Helper method to update the weight of self"""
        if self.weight_type == 'sum':
            self.weight += weight
        else:
            s_w = sum([sub.weight * sub.leafs for sub in self.subtrees])
            self.weight = s_w / self.leafs

    def _leaf_gen(self, value: Any, weight: float) -> None:
        """Creates a leaf node and appends it to self. Correctly updates
        the <self.weight>.
        """
        # Set a flag that detects if the value we are adding is already
        # inside the tree
        found = False
        for subtree in self.subtrees:
            if subtree.value == value:

                # Value presents in Tree
                found = True

                # Update the subtree's weight
                subtree.weight += weight

        # Value does not present in Tree -> Need to add new value
        # *Note this must be added as a leaf.*
        if not found:
            # Create the new SimplePrefixTree
            leaf = CompressedPrefixTree(self.wgt_type)
            leaf.value = value
            leaf.weight = weight
            leaf.leafs = 1

            # Record the current self.leafs for later calculation

            # Append the leaf to self.subtrees
            self.subtrees.append(leaf)

        # Update the attributes of self
        if self.wgt_type == 'sum':
            self.weight += weight
            self.leafs = sum([sub.leafs for sub in self.subtrees])
        else:
            # Critical calculation. Note self.leafs is already updated
            # at this point.
            self.leafs = sum([sub.leafs for sub in self.subtrees])
            s_w = sum([sub.leafs*sub.weight for sub in self.subtrees])
            self.weight = s_w / self.leafs

    def _create_int(self, prefix: List) -> CompressedPrefixTree:
        """Creates a node that will act as an internal node in the tree.
        """
        internal = CompressedPrefixTree(self.weight_type)
        internal.value = prefix

        return internal

    def _find_common2(self, prefix: List) -> List:
        """Returns the longest common prefix of <self.value> and <prefix>.
        """
        if not self.value:
            return []

        common = []
        # if isinstance(self.value, list):
        if not self.is_leaf():
            for i in range(min(len(prefix), len(self.value))):
                if self.value[i] != prefix[i]:
                    break
                common.append(prefix[i])

        return common

    def _find_common(self, prefix: List) -> List:
        """Returns the longest common prefix of <prefix> and all of this tree's
        subtrees.
        """
        common_all = []
        for subtree in self.subtrees:
            # if isinstance(subtree.value, list):
            if not subtree.is_leaf():
                common = []
                for i in range(min(len(prefix), len(subtree.value))):
                    if subtree.value[i] != prefix[i]:
                        break
                    common.append(prefix[i])

                common_all.append(common)

        best_common = []
        for item in common_all:
            if len(item) > len(best_common):
                best_common = item

        return best_common

    def autocomplete(self, prefix: List,
                     limit: Optional[int] = None) -> List[Tuple[Any, float]]:
        """Return up to <limit> matches for the given prefix.

        The return value is a list of tuples (value, weight), and must be
        ordered in non-increasing weight. (You can decide how to break ties.)

        If limit is None, return *every* match for the given prefix.

        Precondition: limit is None or limit > 0.
        """
        values = []
        if prefix == []:
            for subtree in self.subtrees:
                values.extend(subtree._find_leaves())
                if limit is not None and limit <= len(values):
                    break
        else:
            found = False
            # Longest common between self.value and prefix
            self_long_com = self._find_common2(prefix)

            # If the prefix is the self.value, nice, we can directly return
            # everything below
            if prefix == self.value:
                found = True
                values.extend(self.autocomplete([], limit))

            # If the prefix is the same length as long_com, and they are both
            # smaller than self.value, so return everything in self
            elif len(self_long_com) == len(prefix) < len(self.value):
                found = True
                values.extend(self.autocomplete([], limit))

            # Recursion through trees being done here
            else:
                for subtree in self.subtrees:
                    sub_long_com = subtree._find_common2(prefix)
                    if prefix == subtree.value:
                        found = True
                        values.extend(subtree.autocomplete([], limit))
                    elif len(subtree.value) == len(sub_long_com) < len(prefix):
                        found = True
                        values.extend(subtree.autocomplete(prefix, limit))
                    elif len(sub_long_com) == len(prefix) < len(subtree.value):

                        found = True
                        values.extend(subtree.autocomplete([], limit))

            if not found:
                return []

        values = sorted(values,
                        key=lambda sort_wgt: sort_wgt[1], reverse=True)

        if limit is not None:
            return values[0:limit]
        else:
            return values

    def remove(self, prefix: List) -> None:
        """Remove all values that match the given prefix.
        """

        # Base Case, prefix sequence is empty (Self.leafs handled in recursion)
        if prefix == []:
            self.subtrees = []
            self.weight = 0
            self.leafs = 0

        # Prefix sequence non empty -> Recursion
        else:
            # longest no
            long_com = self._find_common2(prefix)

            if long_com == prefix:
                self.remove([])
            else:
                self._remove_helper(prefix)

    def _remove_helper(self, prefix: List) -> None:
        """A helper method to shorten code(?).
        """
        sub_longest_com = self._find_common(prefix)
        for subtree in self.subtrees:
            if subtree._pre_in_sub(sub_longest_com):
                if len(prefix) <= len(subtree.value):
                    subtree.remove([])
                    self.subtrees.remove(subtree)
                elif len(prefix) > len(sub_longest_com) < \
                        len(subtree.value):
                    break
                else:
                    subtree.remove(prefix)

                if len(self.subtrees) == 1:
                    for sub in self.subtrees:
                        self.value = sub.value
                        self.subtrees = sub.subtrees

                self.leafs = sum([s.leafs for s in self.subtrees])

                if self.wgt_type == 'sum':
                    self.weight = sum([s.weight for s in self.subtrees])
                else:
                    self.weight = sum([s.weight * s.leafs for s in
                                       self.subtrees]) / self.leafs
                break


if __name__ == '__main__':
    import python_ta
    python_ta.check_all(config={
        'max-nested-blocks': 4
    })
