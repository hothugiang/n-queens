from dimod import BinaryQuadraticModel

def exact_cover_bqm(problem_set, subsets):
    bqm = BinaryQuadraticModel({}, {}, 0, 'BINARY')

    element_subsets = {}  # Lưu trữ subset mà mỗi phần tử thuộc về

    for i, subset in enumerate(subsets):
        for element in subset:
            if element not in element_subsets:
                element_subsets[element] = set()
            element_subsets[element].add(i)

    for element in problem_set:
        bqm.offset += 1

        if element in element_subsets:
            subsets_containing_element = element_subsets[element]

            for i in subsets_containing_element:
                bqm.add_variable(i, -1)

                for j in range(i):
                    if j in subsets_containing_element:
                        bqm.add_interaction(i, j, 2)

    return bqm

def exact_cover_bqm2(problem_set, subsets):
    bqm = BinaryQuadraticModel({}, {}, 0, 'BINARY')

    element_subsets = {} 

    for i, subset in enumerate(subsets):
        for element in subset:
            if element not in element_subsets:
                element_subsets[element] = set()
            element_subsets[element].add(i)

    for element in problem_set:
        bqm.offset += 1

        if element in element_subsets:
            subsets_containing_element = element_subsets[element]

            for i in subsets_containing_element:
                bqm.add_variable(i, -1)

                for j in range(i):
                    if j in subsets_containing_element:
                        bqm.add_interaction(i, j, 2)

    return bqm