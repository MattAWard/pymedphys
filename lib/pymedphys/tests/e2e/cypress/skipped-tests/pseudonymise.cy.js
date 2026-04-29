/**
 * @license
 * Copyright 2018-2020 Streamlit Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Migrated from: cypress/skipped-tests/pseudonymise.spec.js
// To:            cypress/skipped-tests/pseudonymise.cy.js
//
// Changes for Cypress 15:
//   - Renamed _spec.js → .cy.js (optional for skipped tests, done for consistency)
//   - .attachFile() replaced with .selectFile() (cypress-file-upload removed)
//     .attachFile(fileObj, { subjectType: "drag-n-drop", events: [...] })
//     becomes:
//     .selectFile(fileObj, { action: "drag-drop", force: true })
//

const path = require('path')

const fileName1 = "RS.1.2.840.10008.5.1.4.1.1.481.3.1591744445_Anonymised.dcm";
const fileName2 = "CT.1.3.12.2.1107.5.1.4.115496.30000017121402274359200000404_Anonymised.dcm";

describe("st.file_uploader", () => {
  // with respect to the current working folder
  const downloadsFolder = 'cypress/downloads'

  beforeEach(() => {
    cy.start("pseudonymise")
  });

  it('downloads remote zip', {}, () => {
    // Commented out in original — kept as-is
  })

  it("shows widget correctly", () => {
    cy.get("[data-testid='stFileUploader']")
      .first()
      .should("exist");
    cy.get("[data-testid='stFileUploader'] label")
      .first()
      .should("have.text", "Files to pseudonymise, refresh page after downloading zip(s)");
  });

  it("hides deprecation warning", () => {
    cy.get("[data-testid='stFileUploader']")
      .last()
      .parent()
      .prev()
      .should("not.contain", "FileUploaderEncodingWarning");
  });

  it("uploads single file only", () => {
    // CHANGED: .attachFile() → .selectFile() for Cypress 15
    // cy.selectFile() accepts a path string or a { contents, fileName, mimeType } object.
    // For fixture files, use the path relative to cypress/fixtures/.
    cy.get("[data-testid='stFileUploadDropzone']")
      .eq(0)
      .selectFile(`cypress/fixtures/${fileName1}`, {
        action: "drag-drop",
        force: true,
      });

    cy.get(".uploadedFileName")
      .should("have.text", fileName1);
  });

  it("uploads multiple files", () => {
    // CHANGED: .attachFile() chained calls → single .selectFile() with array
    // cy.selectFile() accepts an array of files for multi-file upload.
    cy.get("[data-testid='stFileUploadDropzone']")
      .eq(0)
      .selectFile(
        [
          `cypress/fixtures/${fileName1}`,
          `cypress/fixtures/${fileName2}`,
        ],
        {
          action: "drag-drop",
          force: true,
        }
      );

    // The widget should show the names of the uploaded files in reverse order
    const filenames = [fileName2, fileName1];
    cy.get(".uploadedFileName").each((uploadedFileName, index) => {
      cy.get(uploadedFileName).should("have.text", filenames[index]);
    });
  });

  it("Pseudonymises data", () => {
    cy.get(".stButton button").contains("Pseudonymise").click({ force: true });
    cy.compute();
  });

});
