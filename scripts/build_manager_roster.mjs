// Build a private Excel request list from the cohort and metadata sources.
import fs from "node:fs/promises";
import { execFileSync } from "node:child_process";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const python = "C:\\Users\\aweyr\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe";
const records = JSON.parse(execFileSync(
  python,
  [path.join(projectRoot, "scripts", "extract_manager_roster.py")],
  { encoding: "utf8", maxBuffer: 4 * 1024 * 1024 },
));

const headers = [
  "Study ID", "Genotype ID", "Brand ID", "Breed", "Sex", "PLINK code",
  "Current outcome", "Template status", "Template elimination reason",
  "Outcome missing", "ID link missing", "Information requested",
  "Updated outcome", "Updated elimination reason", "Confirmed brand ID",
];
const keys = [
  "study_id", "genotype_id", "brand_id", "breed", "sex", "plink_code",
  "current_outcome", "template_status", "template_elimination_reason",
  "outcome_missing", "id_link_missing", "information_requested",
];
const workbook = Workbook.create();

function addRoster(name, rows, tableName) {
  const sheet = workbook.worksheets.add(name);
  sheet.showGridLines = false;
  const matrix = [headers, ...rows.map(row => [...keys.map(key => row[key]), "", "", ""])];
  sheet.getRangeByIndexes(0, 0, matrix.length, headers.length).values = matrix;
  sheet.getRangeByIndexes(0, 0, matrix.length, headers.length).format.font = {
    name: "Arial", size: 10, color: "#243247",
  };
  sheet.getRange("A1:O1").format = {
    fill: "#253D59",
    font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" },
    rowHeight: 30,
    verticalAlignment: "center",
  };
  sheet.getRangeByIndexes(1, 12, rows.length, 3).format.fill = "#FFF2CC";
  sheet.getRangeByIndexes(1, 0, rows.length, headers.length).format.rowHeight = 22;
  const widths = [13, 15, 13, 10, 11, 13, 29, 20, 27, 17, 17, 53, 23, 28, 22];
  widths.forEach((width, index) => {
    sheet.getRangeByIndexes(0, index, matrix.length, 1).format.columnWidth = width;
  });
  sheet.freezePanes.freezeRows(1);
  sheet.tables.add(`A1:O${matrix.length}`, true, tableName);
  return sheet;
}

const requests = records.filter(row => row.information_requested);
addRoster("Request update", requests, "RequestedDogUpdates");
addRoster("All dogs", records, "AllCohortDogs");

workbook.recalculate();
const check = await workbook.inspect({
  kind: "table", range: "Request update!A1:O8", include: "values",
  tableMaxRows: 8, tableMaxCols: 15, maxChars: 6000,
});
console.log(check.ndjson);
const outputDir = path.join(projectRoot, "data", "processed", "private");
await fs.mkdir(outputDir, { recursive: true });
const preview = await workbook.render({
  sheetName: "Request update", range: "A1:O8", scale: 1, format: "png",
});
await fs.writeFile(path.join(outputDir, "manager_roster_preview.png"), new Uint8Array(await preview.arrayBuffer()));
const output = await SpreadsheetFile.exportXlsx(workbook);
const outputPath = path.join(outputDir, "2022TSA_manager_outcome_request.xlsx");
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, allDogs: records.length, requests: requests.length,
  missingOutcomes: records.filter(row => row.outcome_missing === "Yes").length,
  missingLinks: records.filter(row => row.id_link_missing === "Yes").length }));
