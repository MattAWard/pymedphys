// cypress/support/commands.js
//
// Migrated from Cypress 9 to Cypress 15.
//
// CHANGES FROM ORIGINAL:
//   1. Removed `import 'cypress-file-upload'` — replaced by built-in cy.selectFile()
//      in Cypress 12+. The pseudonymise spec needs updating to use cy.selectFile().
//   2. Removed `Cypress.Cookies.defaults({ preserve: ... })` from cy.start() — this API
//      was removed in Cypress 12. In Cypress 12+, cookies are preserved between tests
//      within the same spec by default, so the _xsrf cookie preservation is automatic.
//   3. Updated cy.compute() — `.StatusWidget-enter-done` was the CSS class in older
//      Streamlit. In Streamlit >=1.44, the running indicator may use a different class
//      or data-testid. Added a comment noting this may need updating after Streamlit upgrade.
//   4. All Cypress.Commands.add() calls are unchanged — this API is stable across versions.
//

// import 'cypress-file-upload';
// ^^^ REMOVED: Use cy.selectFile() instead (built-in from Cypress 12+).
//     For the pseudonymise spec, change:
//       .attachFile(file, { force: true, subjectType: "drag-n-drop", events: ["dragenter", "drop"] })
//     to:
//       .selectFile(file, { action: "drag-drop", force: true })

function getBaseUrl() {
  let url = Cypress.env('PYMEDPHYS_GUI_URL')
  if (url === undefined) {
    url = "http://localhost:8501"
  }

  return url
}

Cypress.Commands.add('compute', () => {
  // NOTE: `.StatusWidget-enter-done` is the CSS class used by Streamlit ~1.34.
  // After upgrading Streamlit to >=1.44, inspect the DOM to confirm this class
  // still exists. It may have changed to `[data-testid="stStatusWidget"]` or
  // similar. If so, update the selectors below.
  let start = new Date().getTime();
  cy.get(".StatusWidget-enter-done", { timeout: 4000 }).should($el => {
    let now = new Date().getTime();
    if (now - start < 1000) {
      expect($el).to.exist
    } else {
      expect($el).to.not.exist
    }
  })

  cy.get(".StatusWidget-enter-done", { timeout: 120000 }).should("not.exist", { timeout: 120000 })
})

Cypress.Commands.add('textMatch', (label, length, result) => {
  // From https://github.com/streamlit/streamlit/blob/a03d3b9/e2e/specs/st_markdown.spec.js#L24
  let text = cy.get(`.element-container .stMarkdown p:contains(${label})`).should("have.length", length)
  if (result !== null) {
    text.find('code').each(($el) => {
      return cy.wrap($el).should("have.text", result)
    })
  }
})

Cypress.Commands.add('start', (app) => {
  let url = getBaseUrl()
  cy.visit(url)

  // REMOVED: Cypress.Cookies.defaults({ preserve: ["_xsrf"] })
  // This API was removed in Cypress 12. Cookies are now preserved between
  // tests within the same spec by default in Cypress 12+.

  cy.visit(`${url}/?app=${app}`);
  cy.compute()

  // From https://github.com/streamlit/streamlit/blob/a03d3b9/e2e/specs/component_template.spec.js#L41-L42
  // Make the ribbon decoration line disappear
  cy.get("[data-testid='stDecoration']").invoke("css", "display", "none");

  cy.compute()
})

Cypress.Commands.add('scroll', () => {
  cy.get(".main").scrollTo("bottomLeft");
})

Cypress.Commands.add('finalScreenshot', () => {
  cy.scroll()
  cy.compute()
  cy.scroll()
  cy.compute()
  cy.scroll()
  cy.screenshot()
})

Cypress.Commands.add('radio', (title, item) => {
  cy
    .contains(title)
    .parent()
    .contains(item)
    .find("input")
    .first()
    .click({ force: true });
  cy.compute()
})
