#!/usr/bin/env python
import json
import matplotlib.pyplot as plt

def find_sample_percentages(counts: dict, total):
    percentages = {}
    hopped = 0
    for pair in counts:
        split = pair.split("_")
        idx_1 = split[0]
        idx_2 = split[1]
        
        if idx_1 == idx_2:
            percentages[pair] = (counts[pair] / total) * 100
        else:
            hopped += int(counts[pair])
    return percentages, hopped, hopped/total


def show_stats(percentages, hopped, hopped_perc):
    stats = """"""
    for perc in percentages:
        stats += f"{perc}: {percentages[perc]}\n"
    
    print("The following barcodes matched with the following percentatges")
    print(stats)
    print("--------------------------------------------------------------")
    print()
    
    print(f"The number of indexes that hopped is {hopped} which is {hopped_perc*100}% of all reads")
    
    
    
def plot_percentages(percentages):
    
    x = percentages.keys()
    x = [pair.split("_")[0] for pair in x]
    y = percentages.values()
    for value in y:
        value = float(value)
        
    plt.bar(x, y)
    plt.xticks(rotation=45,ha='right')
    plt.xlabel('Barcodes')
    plt.ylabel('Percent match')
    plt.title('Percent match for each of the sample barcodes')
    plt.tight_layout()
    
    plt.savefig("percentages.png")
    
if __name__ == "__main__":
           
    with open('counts.json', 'r') as fh:
        data = json.load(fh)
    

    total = sum(data.values())
        
    stats, hopped, hopped_perc = find_sample_percentages(data, total)

    show_stats(stats, hopped, hopped_perc) 
    plot_percentages(stats)
