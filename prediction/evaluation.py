import math
import xlwt
from datetime import timedelta

class Evaluator:

    def __init__(self, input_path):
        self.evaluation_results = []
        self.input_path = input_path

    def evaluate(self, emdb_id, predicted_file, gt_file, execution_time):
        """This method finds the closest pred_ca/gt_ca pair in the entire
        structure. Then removes them from the set and continues to find the next
        closest pair until all pairs with a distance of less than 3A have been
        removed."""
        gt_pdb = open(gt_file, 'r')
        gt_ca_atoms = list()
        for line in gt_pdb:
            if line.startswith("ATOM") and line[13:16] == "CA ":
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                gt_ca_atoms.append(tuple([x, y, z]))
        gt_pdb.close()
        native_ca_atoms = len(gt_ca_atoms)
        prediction_file = open(predicted_file, 'r')
        modeled_ca = 0
        pred_ca_atoms = list()
        previous_index = -2
        for line in prediction_file:
            if line.startswith("ATOM") and line[13:16] == "CA ":
                modeled_ca += 1
                index = int(line[23:26])
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                if index != previous_index + 1:
                    pred_ca_atoms.append(list())
                pred_ca_atoms[len(pred_ca_atoms) - 1].append(tuple([x, y, z]))
                previous_index = index
        prediction_file.close()

        # Find the number incorrect
        incorrect = 0
        for partial_set in pred_ca_atoms:
            for pred_ca in partial_set:
                min_dist = 1000000
                for gt_ca in gt_ca_atoms:
                    dist = distance(pred_ca[2], gt_ca[2], pred_ca[1], gt_ca[1], pred_ca[0], gt_ca[0])
                    if dist < min_dist:
                        min_dist = dist
                if min_dist > 3:
                    incorrect += 1

        total_ca = 0
        squared_sum = 0
        for partial_set in pred_ca_atoms:
            # Do the first direction now
            removed_gt_atoms_one = list()
            one_squared_sum = 0
            one_total_ca = 0
            for pred_ca in partial_set:
                min_dist = 1000000
                closest_gt_ca = None
                for gt_ca in gt_ca_atoms:
                    dist = distance(pred_ca[2], gt_ca[2], pred_ca[1], gt_ca[1], pred_ca[0], gt_ca[0])
                    if dist < min_dist:
                        min_dist = dist
                        closest_gt_ca = gt_ca
                if min_dist < 3:
                    gt_ca_atoms.remove(closest_gt_ca)
                    removed_gt_atoms_one.append(closest_gt_ca)
                    one_total_ca += 1
                    one_squared_sum += min_dist ** 2
            # Restore the gt_list
            for removed_ca in removed_gt_atoms_one:
                gt_ca_atoms.append(removed_ca)
            # Now do the other direction
            partial_set.reverse()
            removed_gt_atoms_two = list()
            two_squared_sum = 0
            two_total_ca = 0
            for pred_ca in partial_set:
                min_dist = 1000000
                closest_gt_ca = None
                for gt_ca in gt_ca_atoms:
                    dist = distance(pred_ca[2], gt_ca[2], pred_ca[1], gt_ca[1], pred_ca[0], gt_ca[0])
                    if dist < min_dist:
                        min_dist = dist
                        closest_gt_ca = gt_ca
                if min_dist < 3:
                    gt_ca_atoms.remove(closest_gt_ca)
                    removed_gt_atoms_two.append(closest_gt_ca)
                    two_total_ca += 1
                    two_squared_sum += min_dist ** 2
            # Restore the gt_list
            for removed_ca in removed_gt_atoms_two:
                gt_ca_atoms.append(removed_ca)
            # Now use the better fit
            one_fit = 0 if one_total_ca == 0 else math.sqrt(one_squared_sum / one_total_ca)
            two_fit = 0 if two_total_ca == 0 else math.sqrt(two_squared_sum / two_total_ca)
            if two_total_ca > one_total_ca or (two_total_ca == one_total_ca and two_fit < one_fit):
                for ca in removed_gt_atoms_two:
                    gt_ca_atoms.remove(ca)
                total_ca += two_total_ca
                squared_sum += two_squared_sum
            else:
                for ca in removed_gt_atoms_one:
                    gt_ca_atoms.remove(ca)
                total_ca += one_total_ca
                squared_sum += one_squared_sum

        self.evaluation_results.append(EvaluationResult(emdb_id,
                                                        modeled_ca,
                                                        native_ca_atoms,
                                                        total_ca,
                                                        total_ca / native_ca_atoms,
                                                        math.sqrt(squared_sum / total_ca),
                                                        incorrect,
                                                        execution_time))

    def create_report(self, output_path, execution_time):
        """Creates excel document containing evaluation reports"""
        # Don't create report if there are no evaluation results
        if not self.evaluation_results:
            return

        book = xlwt.Workbook()
        sh = book.add_sheet('results')

        sh.write(0, 0, 'EMDB ID')
        sh.write(0, 1, '# Modeled Ca Atoms')
        sh.write(0, 2, '# Native Ca Atoms')
        sh.write(0, 3, '# Matching Ca Atoms')
        sh.write(0, 4, 'Matching Percentage')
        sh.write(0, 5, 'RMSD')
        sh.write(0, 6, 'Incorrect')
        sh.write(0, 7, 'Execution Time')

        for i in range(len(self.evaluation_results)):
            sh.write(1 + i, 0, self.evaluation_results[i].name)
            sh.write(1 + i, 1, self.evaluation_results[i].num_modeled_ca)
            sh.write(1 + i, 2, self.evaluation_results[i].num_native_ca)
            sh.write(1 + i, 3, self.evaluation_results[i].num_matching_ca)
            sh.write(1 + i, 4, self.evaluation_results[i].matching_ca_per)
            sh.write(1 + i, 5, self.evaluation_results[i].rmsd)
            sh.write(1 + i, 6, self.evaluation_results[i].num_incorrect)
            sh.write(1 + i, 7, str(timedelta(seconds=int(self.evaluation_results[i].execution_time))))

        rmsd_avg = sum(r.rmsd for r in self.evaluation_results) / len(self.evaluation_results)
        matching_ca_per_avg = sum(r.matching_ca_per for r in self.evaluation_results) / len(self.evaluation_results)
        execution_time_avg = sum(r.execution_time for r in self.evaluation_results) / len(self.evaluation_results)
        sh.write(len(self.evaluation_results) + 1, 0, 'Avg.')
        sh.write(len(self.evaluation_results) + 1, 4, matching_ca_per_avg)
        sh.write(len(self.evaluation_results) + 1, 5, rmsd_avg)
        sh.write(len(self.evaluation_results) + 1, 7, str(timedelta(seconds=int(execution_time_avg))))
        sh.write(len(self.evaluation_results) + 2, 0, 'Total')
        sh.write(len(self.evaluation_results) + 2, 7, str(timedelta(seconds=int(execution_time))))

        book.save(output_path + 'results.xls')


class EvaluationResult:

    def __init__(self, name, num_modeled_ca, num_native_ca, num_matching_ca, matching_ca_per, rmsd, num_incorrect,
                 execution_time):
        self.name = name
        self.num_modeled_ca = num_modeled_ca
        self.num_native_ca = num_native_ca
        self.num_matching_ca = num_matching_ca
        self.matching_ca_per = matching_ca_per
        self.rmsd = rmsd
        self.num_incorrect = num_incorrect
        self.execution_time = execution_time


def distance(z1, z2, y1, y2, x1, x2):
    """Calculates Euclidean distance between two points"""
    z_diff = z1 - z2
    y_diff = y1 - y2
    x_diff = x1 - x2
    sum_squares = math.pow(z_diff, 2) + math.pow(y_diff, 2) + math.pow(x_diff, 2)
    return math.sqrt(sum_squares)
